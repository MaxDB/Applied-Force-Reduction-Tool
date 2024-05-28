function [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] =  parse_h_error_data(Static_Data,L_eigenvectors)
num_dofs = Static_Data.Model.num_dof;

r_eigenvectors = Static_Data.Model.reduced_eigenvectors;
h_eigenvectors = [r_eigenvectors,L_eigenvectors]; 
num_h_modes = size(h_eigenvectors,2);

K_array = Static_Data.tangent_stiffness;
is_fe_system = ~isnumeric(K_array);
num_loadcases = size(K_array,3);
mass = Static_Data.Model.mass;
mass_product = mass*h_eigenvectors;

I_h = eye(num_h_modes);

h_coupling_gradient = zeros(num_dofs,num_h_modes,num_loadcases);
h_stiffness = zeros(num_h_modes,num_h_modes,num_loadcases);
for iLoad = 1:num_loadcases
    if is_fe_system
        K_i = K_array.get_matrix(iLoad);
    else
        K_i = K_array(:,:,iLoad);
    end

    perturbation = K_i\mass_product;
    h_i = mass_product'*perturbation;
    theta_hat_i = perturbation - h_eigenvectors*h_i;

    h_coupling_gradient(:,:,iLoad) = theta_hat_i/h_i;
    h_stiffness(:,:,iLoad) =  I_h/h_i;
end

stiffness = Static_Data.Model.stiffness;
perturbation_0 = stiffness\mass_product;
h_0 = mass_product'*perturbation_0;
theta_0 = perturbation_0 - h_eigenvectors*h_0;

h_coupling_gradient_0 = theta_0/h_0;
h_stiffness_0 =  I_h/h_0;
end