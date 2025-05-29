function Static_Data = apply_small_force(Static_Data,validation_modes)


Model = Static_Data.Model;
r_modes = Model.reduced_modes;
num_r_modes = size(r_modes,2);

r_evec = Model.reduced_eigenvectors;
L_evec = Model.low_frequency_eigenvectors;
r_evec = load_data(r_evec);
L_evec = load_data(L_evec);


all_L_modes = Model.low_frequency_modes;
[~,L_map] = ismember(validation_modes,all_L_modes);
L_evec = L_evec(:,L_map);

h_modes = [r_modes,validation_modes];
num_h_modes = size(h_modes,2);
h_evec = [r_evec,L_evec];

stiffness_pointer = Static_Data.get_dataset_values("tangent_stiffness");
is_fe_system = ~isnumeric(stiffness_pointer);
num_loadcases = size(stiffness_pointer,3);
mass = Model.mass;
h_disp_transform = h_evec'*mass;

lambda_all = select_perturbation_scale_factor(Model);
h_map = [1:num_r_modes,L_map+num_r_modes];
lambda = lambda_all(h_map);
Static_Data.perturbation_scale_factor = lambda;

F_h = lambda.*eye(num_h_modes);


perturbation_disp = zeros([size(h_disp_transform'),num_loadcases]);

applied_force = h_disp_transform'*F_h;

parfor iLoad = 1:num_loadcases
    if is_fe_system
        K_i = stiffness_pointer.get_matrix(iLoad);
    else
        K_i = stiffness_pointer(:,:,iLoad);
    end
    perturbation_disp(:,:,iLoad) = K_i\applied_force;
end

data_path = get_data_path(Static_Data) + "perturbation";
Static_Data.perturbation_displacement = Perturbation_Pointer(Model,perturbation_disp,"path",data_path);
end


