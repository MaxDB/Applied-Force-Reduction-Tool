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
% num_h_modes = size(L_modes,2);
% h_evec = [L_evec];
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

%%% TEST
% Plot_Settings.plot_type = "physical";
% Plot_Settings.coords = [3,2,1];
% Manifold_One.system = "mass_spring_roller_1";
% Manifold_Two.system = "mass_spring_roller_12";
% ax = compare_stress_manifold({Manifold_One,Manifold_Two},"opts",Plot_Settings);
% hold(ax,"on")
% x_tilde = Static_Data.get_dataset_values("physical_displacement");
%%%%

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


    %%% TEST
    % x_tilde_i = x_tilde(:,iLoad);
    % x_i = disp_hat_i + x_tilde_i;
    % x_minus_i = -disp_hat_i + x_tilde_i;
    % % plot3([x_minus_i(3,1),x_i(3,1)],[x_minus_i(2,1),x_i(2,1)],[x_minus_i(1,1),x_i(1,1)])
    % plot3([x_minus_i(3,2),x_i(3,2)],[x_minus_i(2,2),x_i(2,2)],[x_minus_i(1,2),x_i(1,2)])
end

stiffness = Model.stiffness;
perturbation_0 = stiffness\mass_product;
h_0 = mass_product'*perturbation_0;
disp_0 = perturbation_0;

h_coupling_gradient_0 = disp_0/h_0;
h_stiffness_0 =  I_h/h_0;
end