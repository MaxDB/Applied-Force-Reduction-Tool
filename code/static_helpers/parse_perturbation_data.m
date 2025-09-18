function [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_perturbation_data(Static_Data,L_modes)
CLEAN_DATA = 1;
MIN_DISP = 1e-15;

Model = Static_Data.Model;
r_modes = Model.reduced_modes;
num_r_modes = size(r_modes,2);
r_evec = Model.reduced_eigenvectors;
L_evec = Model.low_frequency_eigenvectors;
r_evec = load_data(r_evec);
L_evec = load_data(L_evec);

all_L_modes = Model.low_frequency_modes;
[~,L_map] = ismember(L_modes,all_L_modes);
L_evec = L_evec(:,L_map);

h_modes = [r_modes,L_modes];
num_h_modes = size(h_modes,2);
h_evec = [r_evec,L_evec];
h_map = [1:num_r_modes,L_map+num_r_modes];


num_dofs = Model.num_dof;
mass = Model.mass;
h_disp_transform = h_evec'*mass;

lambda = Static_Data.perturbation_scale_factor(1,h_map);
F_h = lambda.*eye(num_h_modes);

Perturbation = Static_Data.get_dataset_values("perturbation_displacement");
perturbation_disp = Perturbation.get_displacement(h_modes);


num_loadcases = size(perturbation_disp,3);


h_coupling_gradient = zeros(num_dofs,num_h_modes,num_loadcases);
h_stiffness = zeros(num_h_modes,num_h_modes,num_loadcases);

h_Disp_Transform_Const = parallel.pool.Constant(h_disp_transform);
Perturbation_Disp_Const = parallel.pool.Constant(perturbation_disp);
parfor (iLoad = 1:num_loadcases,get_current_parallel_jobs)
    % for iLoad = 1:num_loadcases

    disp_hat = Perturbation_Disp_Const.Value(:,:,iLoad);
    if CLEAN_DATA
        disp_hat(mean(abs(disp_hat),2) < MIN_DISP,:) = 0;
    end
    h_disp = h_Disp_Transform_Const.Value*disp_hat;

    h_coupling_gradient(:,:,iLoad) = disp_hat/h_disp;

    h_stiffness(:,:,iLoad) =  F_h/h_disp;

end

stiffness = Model.stiffness;
disp_hat_0 = lambda.*(stiffness\(h_disp_transform'));
h_0 = h_disp_transform*disp_hat_0;

h_coupling_gradient_0 = disp_hat_0/h_0;
h_stiffness_0 =  F_h/h_0;
end