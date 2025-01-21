function [orbit_stability,orbit_evals] = get_h_stability(h_terms,t0)
[h_inertia,h_conv,h_stiff,~] = h_terms{:};
num_h_modes = size(h_inertia,1);
num_time_points = size(t0,2);
I_L = eye(num_h_modes);

h_disp_span = 1:num_h_modes;
h_vel_span = h_disp_span + num_h_modes;

orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes,num_time_points);
for iTime = 1:(num_time_points)
    % orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes);
    orbit_jacobian(h_disp_span,h_vel_span,iTime) = I_L;
    orbit_jacobian(h_vel_span,h_disp_span,iTime) = -h_inertia(:,:,iTime)\h_stiff(:,:,iTime);
    orbit_jacobian(h_vel_span,h_vel_span,iTime) = -h_inertia(:,:,iTime)\h_conv(:,:,iTime);
end
fundamental_eq = @(t,z) orbit_jacobian_func(t,z,t0,orbit_jacobian);
[~,fundamental_mat] = ode45(fundamental_eq,[0,t0(end)],eye(2*num_h_modes));
% 
% figure;
% plot(t_ode,fundamental_mat(:,1))

monodromy_mat = reshape(fundamental_mat(end,:)',2*num_h_modes,2*num_h_modes);
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