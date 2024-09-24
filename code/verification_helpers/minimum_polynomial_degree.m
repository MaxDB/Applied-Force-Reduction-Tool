function [degree,max_error_index] = minimum_polynomial_degree(Static_Data,poly_name,input_data,output_data,degree_range)
minimum_degree_start = tic;

max_degree = max_polynomial_degree(poly_name,degree_range(2),input_data);
max_fitting_error = Static_Data.Validation_Options.maximum_fitting_error;

degree = degree_range(1);
rom_degree = ones(4,1);

poly_name = convertStringsToChars(poly_name);
poly_name = string([upper(poly_name(1)),poly_name(2:end)]);
switch poly_name
    case "Force"
        degree_index = 1;
    case "Physical_Displacement"
        degree_index = 2;
        min_disp = Reduced_System(Static_Data,[1,1]).minimum_displacement;
        max_output = max(abs(output_data),[],2);
        output_data(max_output < min_disp,:) = 0;
end
while degree <= max_degree
    rom_degree(degree_index) = degree; 
    Rom = Reduced_System(Static_Data,rom_degree);
    Poly = Rom.(poly_name + "_Polynomial");
    poly_data = Poly.evaluate_polynomial(input_data);
    poly_error = coeff_of_determination(output_data,poly_data);
    [max_poly_error,max_error_index] = max(poly_error);

    if max_poly_error < max_fitting_error
        break
    end

    degree = degree + 2;
end
if degree > max_degree
    degree = max_degree;
    warning(poly_name + " polynomial requires a degree greater than maximum: " + max_degree + ". Worst index is " + max_error_index)
end

minimum_degree_time = toc(minimum_degree_start);
log_message = sprintf(poly_name + " degree is %i: %.1f seconds" ,[degree,minimum_degree_time]);
logger(log_message,3)

end

function max_degree = max_polynomial_degree(poly_name,constrained_max_degree,input_data)
num_points = size(input_data,2);
num_modes = size(input_data,1);
switch poly_name
    case "Force"
        degree_shift = 1;
    case "Physical_Displacement"
        degree_shift = 0;
end

num_coeffs = inf;
poly_order_range = @(max_degree) [3,max_degree + degree_shift + 2];
max_degree = constrained_max_degree + 2;
while num_coeffs > num_points
    max_degree = max_degree - 2;
    num_coeffs = get_num_coeffs(poly_order_range(max_degree),num_modes);
end

end

function num_coeffs = get_num_coeffs(poly_order_range,num_inputs)
    num_coeffs = 0;
    for iDegree = poly_order_range(1):poly_order_range(2)
        num_coeffs = num_coeffs + nchoosek(iDegree + num_inputs - 1,num_inputs - 1);
    end
end