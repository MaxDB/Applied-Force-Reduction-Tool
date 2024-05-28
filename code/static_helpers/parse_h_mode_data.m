function [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_h_mode_data(Static_Data,h_eigenvectors)
num_dofs = Static_Data.Model.num_dof;
num_h_modes = size(h_eigenvectors,2);

K_array = Static_Data.tangent_stiffness;
is_fe_system = ~isnumeric(K_array);
num_loadcases = size(K_array,3);
mass = Static_Data.Model.mass;
mass_product = mass*h_eigenvectors;
I_L = eye(num_h_modes);

r = Static_Data.reduced_displacement;
% num_r_modes = size(r,1);
r_eigenvectors = Static_Data.Model.reduced_eigenvectors;
r_disp_transform = r_eigenvectors'*mass;

theta = Static_Data.condensed_displacement;

% total_loadcases = num_h_modes*num_loadcases;
% r_hat = zeros(num_r_modes,total_loadcases);
% h_hat = zeros(num_h_modes,total_loadcases);
% theta_hat = zeros(num_dofs,total_loadcases);
% f_h_hat = zeros(num_h_modes,total_loadcases);

h_coupling_gradient = zeros(num_dofs,num_h_modes,num_loadcases);
h_stiffness = zeros(num_h_modes,num_h_modes,num_loadcases);

for iLoad = 1:num_loadcases
    if is_fe_system
        K_i = K_array.get_matrix(iLoad);
    else
        K_i = K_array(:,:,iLoad);
    end

    % load_span = ((iLoad-1)*num_h_modes + 1):(num_h_modes*iLoad);

    perturbation = K_i\mass_product + r_eigenvectors*r(:,iLoad) + theta(:,iLoad);
    h_hat_i = mass_product'*perturbation;
    r_hat_i = r_disp_transform*perturbation;

    theta_hat_i = perturbation - h_eigenvectors*h_hat_i - r_eigenvectors*r_hat_i;

    % h_hat(:,load_span) = h_hat_i;
    % r_hat(:,load_span) = r_hat_i;
    % theta_hat(:,load_span) = theta_hat_i;
    % 
    % f_h_hat(:,load_span) = I_L;

    h_coupling_gradient(:,:,iLoad) = theta_hat_i/h_hat_i;
    h_stiffness(:,:,iLoad) = I_L/h_hat_i;
end

stiffness = Static_Data.Model.stiffness;
perturbation_0 = stiffness\mass_product;
h_0 = mass_product'*perturbation_0;
h_coupling_gradient_0 = perturbation_0/h_0 - h_eigenvectors;
h_stiffness_0 =  I_L/h_0;

end