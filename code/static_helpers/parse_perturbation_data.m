function [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_perturbation_data(Static_Data,L_modes)
CLEAN_DATA = 1;
MIN_DISP = 1e-15;

Model = Static_Data.Model;
r_modes = Model.reduced_modes;
num_r_modes = Static_Data.get_reduced_dimension;
r_evec = Model.reduced_eigenvectors;
L_evec = Model.low_frequency_eigenvectors;
r_evec = load_data(r_evec);
L_evec = load_data(L_evec);

all_L_modes = Model.low_frequency_modes;
[~,L_map] = ismember(L_modes,all_L_modes);
L_evec = L_evec(:,L_map);

h_modes = [r_modes,L_modes];
h_evec = [r_evec,L_evec];
h_map = [1:num_r_modes,L_map+num_r_modes];


mapped_h_modes = h_modes(h_modes>2000) - 2000*floor(h_modes(h_modes>2000)/2000);
remove_modes = ismember(h_modes,mapped_h_modes);


h_modes(remove_modes) = [];
h_evec(:,remove_modes) = [];
h_map(remove_modes) = [];

num_h_modes = size(h_modes,2);


num_dofs = Model.num_dof;
mass = Model.mass;
h_disp_transform = h_evec'*mass;

lambda = Static_Data.perturbation_scale_factor(1,h_map);
F_h = lambda.*eye(num_h_modes);

Perturbation = Static_Data.get_dataset_values("perturbation_displacement");
perturbation_disp = Perturbation.get_displacement(h_modes);


num_loadcases = size(perturbation_disp,3);


h_coupling_gradient = zeros(num_dofs,num_h_modes,num_loadcases);
if all(h_modes < 2000)
    h_stiffness = zeros(num_h_modes,num_h_modes,num_loadcases);
else
    h_stiffness = zeros(num_dofs,num_h_modes,num_loadcases);
end


if all(h_modes < 2000)
    h_Disp_Transform_Const = parallel.pool.Constant(h_disp_transform);
    Perturbation_Disp_Const = parallel.pool.Constant(perturbation_disp);

    
    
    parfor (iLoad = 1:num_loadcases,get_current_parallel_jobs)
    % warning("parallelisation disabled")
    % for iLoad = 1:num_loadcases

        disp_hat = Perturbation_Disp_Const.Value(:,:,iLoad);
        if CLEAN_DATA
            disp_hat(mean(abs(disp_hat),2) < MIN_DISP,:) = 0;
        end
        h_disp = h_Disp_Transform_Const.Value*disp_hat;

        h_coupling_gradient(:,:,iLoad) = disp_hat/h_disp;

        h_stiffness(:,:,iLoad) =  F_h/h_disp;

    end
    clear("h_Disp_Transform_Const","Perturbation_Disp_Const")
else
    geometry = load_geometry(Model);
    all_dofs = Model.num_dof + size(Model.dof_boundary_conditions,1);
    node_map = 1:all_dofs;
    node_map(Model.dof_boundary_conditions) = [];


    Mesh_Data = get_mesh_data(geometry);
    Mesh_Data = Mesh_Data{1};

    num_nodes = all_dofs/Mesh_Data.dimension;
    Mesh_Data.num_nodes = num_nodes;

    base_displacements_bc = Static_Data.get_dataset_values("physical_displacement");
    starting_disp = zeros(all_dofs,1);
    base_displacements = zeros(all_dofs,num_loadcases);
    base_displacements(node_map,:) = base_displacements_bc;
    
    applied_force = h_disp_transform'*F_h;

    applied_forces = zeros(all_dofs,num_h_modes);
    applied_forces(node_map,:) = applied_force;
    

    %-----
    % DEBUG
    % figure
    % ax = axes;
    % hold(ax,"on")
    % r = Static_Data.reduced_displacement;
    % 
    % stiffness = Model.stiffness;
    % disp_hat_0 = lambda.*(stiffness\(h_disp_transform'));
    % h_0 = h_disp_transform*disp_hat_0;
    % plot(ax,0,h_0,"rx");
    %-----

    for iLoad = 1:num_loadcases
        
        disp_hat = perturbation_disp(:,:,iLoad);

        base_displacement = base_displacements(:,iLoad);
        rotated_shapes = zeros(size(applied_force));
        for iMode = 1:num_h_modes
            base_shape = applied_forces(:,iMode);
            if h_modes(iMode) > 2000
                rotated_shape = rotate_force(base_shape,starting_disp,base_displacement,Mesh_Data);
            else
                rotated_shape = base_shape;
            end
            rotated_shapes(:,iMode) = rotated_shape(node_map);
        end

        rotated_mode_shape = rotated_shapes./lambda;

        h_disp = rotated_mode_shape'*disp_hat;
        %parametrisation is still dodgy

        h_coupling_gradient(:,:,iLoad) = disp_hat/h_disp;

        h_stiffness(:,:,iLoad) =  applied_force/h_disp;
        %must be incorrect?


        % %-----
        % % DEBUG
        % plot(ax,r(iLoad),0,"kx")
        % plot(ax,r(iLoad),h_disp,"rx")
        % %-----
    end
end


stiffness = Model.stiffness;
disp_hat_0 = lambda.*(stiffness\(h_disp_transform'));
h_0 = h_disp_transform*disp_hat_0;

h_coupling_gradient_0 = disp_hat_0/h_0;

if all(h_modes < 2000)
    h_stiffness_0 =  F_h/h_0;
else
    h_stiffness_0 = applied_force/h_0;
end
end