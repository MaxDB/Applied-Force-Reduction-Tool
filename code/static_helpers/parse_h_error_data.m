function [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] =  parse_h_error_data(Static_Data,L_modes)

Model = Static_Data.Model;
r_modes = Model.reduced_modes;
num_r_modes = size(r_modes,2);
r_evec = Model.reduced_eigenvectors;

all_L_modes = Model.low_frequency_modes;
[~,L_map] = ismember(L_modes,all_L_modes);
L_evec = Model.low_frequency_eigenvectors(:,L_map);

h_modes = [r_modes,L_modes];
num_h_modes = size(h_modes,2);
h_evec = [r_evec,L_evec];
% h_map = [1:num_r_modes,L_map+num_r_modes];

num_dofs = Model.num_dof;

K_array = Static_Data.get_dataset_values("tangent_stiffness");
is_fe_system = ~isnumeric(K_array);
num_loadcases = size(K_array,3);
mass = Model.mass;
mass_product = mass*h_evec;

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
    disp_hat_i = perturbation;

    h_coupling_gradient(:,:,iLoad) = disp_hat_i/h_i;
    h_stiffness(:,:,iLoad) =  I_h/h_i;
end

stiffness = Model.stiffness;
perturbation_0 = stiffness\mass_product;
h_0 = mass_product'*perturbation_0;
disp_0 = perturbation_0;

h_coupling_gradient_0 = disp_0/h_0;
h_stiffness_0 =  I_h/h_0;
end