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
h_disp_transform = h_evec'*mass;

lambda_all = select_perturbation_scale_factor(Model);
h_map = [1:num_r_modes,L_map+num_r_modes];
lambda = lambda_all(h_map);

F_h = lambda.*eye(num_h_modes);

h_coupling_gradient = zeros(num_dofs,num_h_modes,num_loadcases);
h_stiffness = zeros(num_h_modes,num_h_modes,num_loadcases);
%%%test
% perturbation = zeros([size(h_disp_transform'),num_loadcases]);
%%%
for iLoad = 1:num_loadcases
    if is_fe_system
        K_i = K_array.get_matrix(iLoad);
    else
        K_i = K_array(:,:,iLoad);
    end
    disp_hat = K_i\(h_disp_transform'*F_h);
    %%%test
    % perturbation(:,:,iLoad) = disp_hat;
    %%%

    h_disp = h_disp_transform*disp_hat;

    h_coupling_gradient(:,:,iLoad) = disp_hat/h_disp;
    h_stiffness(:,:,iLoad) =  F_h/h_disp;
end

%%%test
% Static_Data.perturbation_displacement = perturbation;
% Static_Data.perturbation_scale_factor = select_perturbation_scale_factor(Model);
% 
% Static_Data_P = load_static_data("clamped_beam_copy_1");
% 
% plot_outputs = [randi(1434,[9,1]),randi(10,[9,1])];
% 
% ax = plot_static_data("perturbation",Static_Data,"outputs",plot_outputs);
% plot_static_data("perturbation",Static_Data_P,"outputs",plot_outputs,"axes",ax);
%%%


stiffness = Model.stiffness;
disp_hat_0 = stiffness\(h_disp_transform'*F_h);
h_disp_0 = h_disp_transform*disp_hat_0;

h_coupling_gradient_0 = disp_hat_0/h_disp_0;
h_stiffness_0 =  F_h/h_disp_0;
end