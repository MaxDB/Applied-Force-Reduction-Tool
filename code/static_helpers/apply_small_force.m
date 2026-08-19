function Static_Data = apply_small_force(Static_Data,validation_modes)

Model = Static_Data.Model;
r_modes = Model.reduced_modes;
num_r_modes = Static_Data.get_reduced_dimension;

r_evec = Model.reduced_eigenvectors;
L_evec = Model.low_frequency_eigenvectors;
r_evec = load_data(r_evec);
L_evec = load_data(L_evec);

h_modes = [r_modes,validation_modes];

all_L_modes = Model.low_frequency_modes;
all_h_modes = [r_modes,all_L_modes];
[~,h_map] = ismember(h_modes,all_h_modes);

h_evec = [r_evec,L_evec];
h_evec = h_evec(:,h_map);


num_h_modes = size(h_modes,2);


stiffness_pointer = Static_Data.get_dataset_values("tangent_stiffness");
is_fe_system = ~isnumeric(stiffness_pointer);
num_loadcases = size(stiffness_pointer,3);
mass = Model.mass;
h_disp_transform = h_evec'*mass;

lambda_all = select_perturbation_scale_factor(Model);
lambda = lambda_all(h_map);
Static_Data.perturbation_scale_factor = lambda;

F_h = lambda.*eye(num_h_modes);


perturbation_disp = zeros([size(h_disp_transform'),num_loadcases]);

applied_force = h_disp_transform'*F_h;


if all(h_modes < 2000)
    Const_Applied_Force = parallel.pool.Constant(applied_force);
    Const_Stiffness = parallel.pool.Constant(stiffness_pointer);


    parfor (iLoad = 1:num_loadcases,get_current_parallel_jobs)
        if is_fe_system
            perturbation_disp(:,:,iLoad) = Const_Stiffness.Value.get_matrix(iLoad)\Const_Applied_Force.Value;
        else
            perturbation_disp(:,:,iLoad) = Const_Stiffness.Value(:,:,iLoad)\Const_Applied_Force.Value;
        end
    end
    clear("Const_Applied_Force","Const_Stiffness")
else

    
    geometry = load_geometry(Model);
    all_dofs = Model.num_dof + size(Model.dof_boundary_conditions,1);
    node_map = 1:all_dofs;
    node_map(Model.dof_boundary_conditions) = [];


    Mesh_Data = get_mesh_data(geometry);
    Mesh_Data = Mesh_Data{1};

    num_nodes = all_dofs/Mesh_Data.dimension;
    Mesh_Data.num_nodes = num_nodes;
    
    displacements_bc = Static_Data.get_dataset_values("physical_displacement");
    starting_disp = zeros(all_dofs,1);
    displacements = zeros(all_dofs,num_loadcases);
    displacements(node_map,:) = displacements_bc;

    applied_forces = zeros(all_dofs,num_h_modes);
    applied_forces(node_map,:) = applied_force;

    for iLoad = 1:num_loadcases
        if is_fe_system
            stiffness = stiffness_pointer.get_matrix(iLoad);
        else
            stiffness = stiffness_pointer(:,:,iLoad);
        end
        
        displacement = displacements(:,iLoad);
        rotated_shapes = zeros(size(applied_force));
        for iMode = 1:num_h_modes
            base_shape = applied_forces(:,iMode);
            if h_modes(iMode) > 2000
                rotated_shape = rotate_force(base_shape,starting_disp,displacement,Mesh_Data);
            else
                rotated_shape = base_shape;
            end
            rotated_shapes(:,iMode) = rotated_shape(node_map);
        end

        perturbation_disp(:,:,iLoad) = stiffness\rotated_shapes;
    end
end


data_path = get_data_path(Static_Data) + "perturbation";
Static_Data.perturbation_displacement = Perturbation_Pointer(Model,perturbation_disp,"path",data_path);
end


