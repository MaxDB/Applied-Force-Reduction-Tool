function Static_Data = minimum_h_degree(Static_Data)
MIN_VALIDATION_RATIO = 1e-1;

num_r_modes = length(Static_Data.Model.reduced_modes);
switch num_r_modes
    case 1
        max_degree = max(15,max(Static_Data.validated_degree([1,2]))+1);
    otherwise
        max_degree = max(Static_Data.validated_degree([1,2]))+1;
end

Validation_Opts = Static_Data.Validation_Options;
max_fitting_error = Validation_Opts.maximum_fitting_error*10;


V = Static_Data.potential_energy;
validated_points = V < Static_Data.Model.energy_limit;

r = Static_Data.reduced_displacement(:,validated_points);

%find minimum h_stiffness degree
h_stiffness = Static_Data.low_frequency_stiffness(:,:,validated_points);
h_stiffness_max = max(abs(h_stiffness),[],3);
h_stiffness_col_norm = h_stiffness_max./max(h_stiffness_max,[],1);
lin_index = find(h_stiffness_col_norm>MIN_VALIDATION_RATIO);

num_h_modes = size(h_stiffness,1);
[row,col] = ind2sub([num_h_modes,num_h_modes],lin_index);
modelled_output_indicies = [row,col];
num_modelled_outputs = size(modelled_output_indicies,1);

h_stiffness_degree = Static_Data.validated_degree(1);
while h_stiffness_degree  < max_degree
    rom = Reduced_System(Static_Data,[1,1,h_stiffness_degree,1]);
    
    stiffness_error = zeros(num_modelled_outputs,1);
    for iOutput = 1:num_modelled_outputs
        h_stiffness_i = squeeze(h_stiffness(modelled_output_indicies(iOutput,1),modelled_output_indicies(iOutput,2),:))';
        h_stiffness_rom = rom.Low_Frequency_Stiffness_Polynomial.evaluate_polynomial(r,modelled_output_indicies(iOutput,:));
        stiffness_error(iOutput) = coeff_of_determination(h_stiffness_i,h_stiffness_rom);
    end
    [max_stiffness_error,max_index] = max(stiffness_error);
    max_stiffness_index = modelled_output_indicies(max_index,:); 
    if max_stiffness_error < max_fitting_error
        break
    end

    h_stiffness_degree = h_stiffness_degree + 2;
end
if h_stiffness_degree >= max_degree
    warning("Low frequency stiffness polynomial requires too high a degree: worst index is " + ...
        "[" + max_stiffness_index(1) + "," + max_stiffness_index(2) + "]")
end
Static_Data.validated_degree(3) = h_stiffness_degree;


%find minimum h_coupling_gradient degree
h_coupling_gradient = Static_Data.low_frequency_coupling_gradient(:,:,validated_points);
h_coupling_gradient_max = max(abs(h_coupling_gradient),[],3);
h_coupling_gradient_col_norm = h_coupling_gradient_max./max(h_coupling_gradient_max,[],1);
lin_index = find(h_coupling_gradient_col_norm>MIN_VALIDATION_RATIO);

num_dofs = size(h_coupling_gradient,1);
num_points = size(h_coupling_gradient,3);

h_coupling_gradient_vector = reshape(h_coupling_gradient,num_dofs*num_h_modes,num_points);
h_coupling_gradient_i = h_coupling_gradient_vector(lin_index,:);

[row,col] = ind2sub([num_dofs,num_h_modes],lin_index);
modelled_output_indicies = [row,col];

h_coupling_gradient_degree = Static_Data.validated_degree(2);
while h_coupling_gradient_degree  < max_degree
    rom = Reduced_System(Static_Data,[1,1,1,h_coupling_gradient_degree]);
    
    h_coupling_gradient_rom = rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r,lin_index);
    coupling_gradient_error = coeff_of_determination(h_coupling_gradient_i,h_coupling_gradient_rom); 

    [max_coupling_gradient_error,max_index] = max(coupling_gradient_error);
    max_coupling_gradient_index = modelled_output_indicies(max_index,:); 
    if max_coupling_gradient_error < max_fitting_error
        break
    end
    modelled_output_indicies(coupling_gradient_error < max_fitting_error,:) = [];
    h_coupling_gradient_degree = h_coupling_gradient_degree + 2;
end
if h_coupling_gradient_degree >= max_degree
    warning("Low frequency coupling gradient polynomial requires too high a degree: worst index is " + ...
        "[" + max_coupling_gradient_index(1) + "," + max_coupling_gradient_index(2) + "]")
end
Static_Data.validated_degree(4) = h_coupling_gradient_degree;
end