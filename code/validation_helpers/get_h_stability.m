function [orbit_stability,orbit_evals] = get_h_stability(h_terms,t0,num_orbit_harmonics)
% STABILITY_METHOD = "time";
STABILITY_METHOD = "koopman-hill";
[h_inertia,h_conv,h_stiff,~] = h_terms{:};
num_h_modes = size(h_inertia,1);

I_L = eye(num_h_modes);

h_disp_span = 1:num_h_modes;
h_vel_span = h_disp_span + num_h_modes;

switch STABILITY_METHOD
    case "time"
        % num_time_points = size(t0,2);
        % orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes,num_time_points);
        % for iTime = 1:(num_time_points)
        %     % orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes);
        %     orbit_jacobian(h_disp_span,h_vel_span,iTime) = I_L;
        %     orbit_jacobian(h_vel_span,h_disp_span,iTime) = -h_inertia(:,:,iTime)\h_stiff(:,:,iTime);
        %     orbit_jacobian(h_vel_span,h_vel_span,iTime) = -h_inertia(:,:,iTime)\h_conv(:,:,iTime);
        % end
        % ode_opts = odeset("RelTol",1e-9,"AbsTol",1e-11);
        % fundamental_eq = @(t,z) orbit_jacobian_func(t,z,t0,orbit_jacobian);
        % [~,fundamental_mat] = ode45(fundamental_eq,[0,t0(end)],eye(2*num_h_modes),ode_opts);
        % monodromy_mat = reshape(fundamental_mat(end,:)',2*num_h_modes,2*num_h_modes);

        num_harmonics = num_orbit_harmonics;
        num_time_points = 1+8*num_harmonics;
        if num_time_points > length(t0)-1
            num_time_points = length(t0)-1;
            if mod(num_time_points,2) == 0
                num_time_points = num_time_points - 1;
            end
        end

        time_points = linspace(0,t0(end),num_time_points+1);
        time_points(end) = [];


        interp_inertia = zeros(num_h_modes,num_h_modes,num_time_points);
        interp_conv = zeros(num_h_modes,num_h_modes,num_time_points);
        interp_stiff = zeros(num_h_modes,num_h_modes,num_time_points);
        for iRow = 1:num_h_modes
            for iCol = 1:num_h_modes
                interp_inertia(iRow,iCol,:) = interp1(t0,squeeze(h_inertia(iRow,iCol,:)),time_points);
                interp_conv(iRow,iCol,:) = interp1(t0,squeeze(h_conv(iRow,iCol,:)),time_points);
                interp_stiff(iRow,iCol,:) = interp1(t0,squeeze(h_stiff(iRow,iCol,:)),time_points);
            end
        end


        orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes,num_time_points);
        for iTime = 1:num_time_points
            % orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes);
            orbit_jacobian(h_disp_span,h_vel_span,iTime) = I_L;
            orbit_jacobian(h_vel_span,h_disp_span,iTime) = -interp_inertia(:,:,iTime)\interp_stiff(:,:,iTime);
            orbit_jacobian(h_vel_span,h_vel_span,iTime) = -interp_inertia(:,:,iTime)\interp_conv(:,:,iTime);
        end
        fourier_coefficients = fft(orbit_jacobian,size(orbit_jacobian,3),3)/num_time_points;

        
        period = t0(end);
        omega = 2*pi/period;
        fundamental_eq = @(t,z) orbit_jacobian_fourier_func(t,z,omega,fourier_coefficients);
        [~,fundamental_mat] = ode45(fundamental_eq,[0,period],eye(2*num_h_modes));
        monodromy_mat = reshape(fundamental_mat(end,:)',2*num_h_modes,2*num_h_modes);



    case "koopman-hill"

        period = t0(end);

        num_harmonics = num_orbit_harmonics;
        num_time_points = 1+8*num_harmonics;
        if num_time_points > length(t0)-1
            num_time_points = length(t0)-1;
            if mod(num_time_points,2) == 0
                num_time_points = num_time_points - 1;
            end
        end


        % num_time_points = size(t0,2) - 1;
        time_points = linspace(0,t0(end),num_time_points+1);
        time_points(end) = [];


        interp_inertia = zeros(num_h_modes,num_h_modes,num_time_points);
        interp_conv = zeros(num_h_modes,num_h_modes,num_time_points);
        interp_stiff = zeros(num_h_modes,num_h_modes,num_time_points);
        for iRow = 1:num_h_modes
            for iCol = 1:num_h_modes
                interp_inertia(iRow,iCol,:) = interp1(t0,squeeze(h_inertia(iRow,iCol,:)),time_points);
                interp_conv(iRow,iCol,:) = interp1(t0,squeeze(h_conv(iRow,iCol,:)),time_points);
                interp_stiff(iRow,iCol,:) = interp1(t0,squeeze(h_stiff(iRow,iCol,:)),time_points);
            end
        end


        orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes,num_time_points);
        for iTime = 1:num_time_points
            % orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes);
            orbit_jacobian(h_disp_span,h_vel_span,iTime) = I_L;
            orbit_jacobian(h_vel_span,h_disp_span,iTime) = -interp_inertia(:,:,iTime)\interp_stiff(:,:,iTime);
            orbit_jacobian(h_vel_span,h_vel_span,iTime) = -interp_inertia(:,:,iTime)\interp_conv(:,:,iTime);
        end
        fourier_coefficients = fft(orbit_jacobian,size(orbit_jacobian,3),3)/num_time_points;


        shifted_coeffs = fftshift(fourier_coefficients,3);

        hill_matrix = get_hill_matrix(shifted_coeffs,period);


        matrix_size = size(fourier_coefficients,1);
        num_matrices = size(hill_matrix,1)/matrix_size;
        max_harmonics = (num_matrices-1)/2;

        left_lin_map = [zeros(matrix_size,matrix_size*max_harmonics),eye(matrix_size),zeros(matrix_size,matrix_size*max_harmonics)];
        right_lin_map = repmat(eye(matrix_size),[num_matrices,1]);

        monodromy_mat = left_lin_map*expm(hill_matrix*period)*right_lin_map;

end


orbit_evals = eig(monodromy_mat);
orbit_stability = max(abs(orbit_evals));

end

function dz = orbit_jacobian_func(t,z,t0,orbit_jacobian)
num_h_modes = size(orbit_jacobian,1)/2;
h_span = 1:num_h_modes;
h_dot_span = h_span + num_h_modes;

approx_jacobian = zeros(2*num_h_modes);
approx_jacobian(h_span,h_dot_span) = eye(num_h_modes);
for iRow = h_dot_span
    for iCol = 1:(2*num_h_modes)
        approx_jacobian(iRow,iCol) = interp1(t0,squeeze(orbit_jacobian(iRow,iCol,:)),t);
    end
end

dz = approx_jacobian*reshape(z,2*num_h_modes,2*num_h_modes);
dz = reshape(dz,4*num_h_modes^2,1);
end
%------------------$
function dz = orbit_jacobian_fourier_func(t,z,frequency,jacobian_coeffs)
num_h_modes = size(jacobian_coeffs,1)/2;

approx_jacobian = evaluate_fourier(jacobian_coeffs,t,frequency);

dz = approx_jacobian*reshape(z,2*num_h_modes,2*num_h_modes);
dz = reshape(dz,4*num_h_modes^2,1);
end

function X = evaluate_fourier(coeffs,t,frequency)
num_harmonics = (size(coeffs,3) - 1)/2;
zero_point = 1;
X = coeffs(:,:,zero_point);

for k = 1:num_harmonics
    % X = X + coeffs(:,:,mid_point+k)*exp(1i*frequency*k*t) + coeffs(:,:,mid_point-k)*exp(-1i*frequency*k*t);
    X = X + 2*real(coeffs(:,:,zero_point+k))*cos(frequency*k*t) + 2*imag(coeffs(:,:,zero_point+k))*sin(frequency*k*t);
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function hill_matrix = get_hill_matrix(fourier_coefficients,period)

matrix_size = size(fourier_coefficients,1);

num_coeffs = size(fourier_coefficients,3);
num_harmonics = (num_coeffs-1)/2;
max_harmonic = floor(num_harmonics/2);

num_usable_coeffs = 2*max_harmonic + 1;

hill_size = matrix_size*(num_usable_coeffs);
hill_matrix = zeros(hill_size);

frequency = 2*pi/period;
frequency_matrix = frequency*eye(matrix_size);




for iRow = 1:num_usable_coeffs
    row_span = (matrix_size*(iRow-1) + 1):(matrix_size*iRow);
    index = 0+iRow;
    
    for iCol = 1:num_usable_coeffs
        col_span = (matrix_size*(iCol-1) + 1):(matrix_size*iCol);
        index = index - 1;
       
        fourier_matrix = get_fourier_matrix(fourier_coefficients,index);
        % fourier_matrix(:,:) = index;
        % frequency_matrix(:,:) = 1;
        if iRow == iCol
            fourier_matrix = fourier_matrix + 1i*frequency_matrix*(max_harmonic+1 - iRow);
        end
        hill_matrix(row_span,col_span) = fourier_matrix;
    end
end
end
%--------------------------------------------------------------------------
function fourier_matrix = get_fourier_matrix(fourier_coefficients,index)
% fourier coefficients are shifted
num_coefficients = size(fourier_coefficients,3);
mid_point = ceil((num_coefficients+1)/2);
fourier_matrix = fourier_coefficients(:,:,mid_point + index);
end