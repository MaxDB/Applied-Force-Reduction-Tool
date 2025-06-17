function max_fitting_degree = get_max_poly_degree(type,num_r_modes,num_points,max_degree)
degrees = 2:1:max_degree;
num_degrees = size(degrees,2);

if num_points == 0
    max_fitting_degree = 1;
    return
end

switch type
    case "force"
        %must be odd
        degree_shift = 1;
        step_back = @(n_fail) 3 + mod(n_fail,2);
        information_per_point = num_r_modes;
    case "displacement"
        % step_back = @(n_fail) 3;
        step_back = @(n_fail) 3 + mod(n_fail,2);
        degree_shift = 0;
        information_per_point = 1;
end

num_coeffs = 0;
for iTerm = 1:num_degrees
    degree = degrees(iTerm);
    num_coeffs = num_coeffs + nchoosek(degree + degree_shift +(num_r_modes-1),num_r_modes-1);
    required_points = ceil(num_coeffs/information_per_point);
    if required_points > num_points
        max_fitting_degree = degrees(iTerm-step_back(degree));
        break
    elseif iTerm == num_degrees
        max_fitting_degree = degrees(end + 1 - step_back(degree+1));
    end
end

end