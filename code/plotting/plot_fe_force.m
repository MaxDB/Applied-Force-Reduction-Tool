function ax = plot_fe_force(Model,displacement,force,varargin)
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%



geometry = Model.load_geometry;
Mesh_Data = get_mesh_data(geometry);
Mesh_Data = Mesh_Data{1};


if Mesh_Data.type ~= "beam"
    return
end

if Mesh_Data.model_dimension ~= 3
    return
end
disp_dofs = 1:3;
rotation_dofs = 4:6;


base_node_position = read_abaqus_node_position(geometry)';
element_members = read_abaqus_element_membership(geometry);


num_nodes = size(base_node_position,2);

all_dofs = Model.num_dof + numel(Model.dof_boundary_conditions);
node_map = 1:all_dofs;
node_map(Model.dof_boundary_conditions) = [];

if size(displacement,1) == size(node_map,2)
    displacement_bc = displacement;
    displacement = zeros(all_dofs,1);
    displacement(node_map) = displacement_bc;
end

if size(force,1) == size(node_map,2)
    force_bc = force;
    force = zeros(all_dofs,1);
    force(node_map) = force_bc;
end

nodal_displacement = reshape(displacement',Mesh_Data.dimension,num_nodes);
nodal_position = base_node_position + nodal_displacement(disp_dofs,:);
nodal_rotation = nodal_displacement(rotation_dofs,:);

nodal_force_all = reshape(force',Mesh_Data.dimension,num_nodes);
nodal_force = nodal_force_all(disp_dofs,:);
nodal_moments = nodal_force_all(rotation_dofs,:);

moment_factor = 1e6;
nodal_force_end = nodal_position + nodal_force;
nodal_moment_end = nodal_position + moment_factor*nodal_moments;


if isempty(ax)
    figure
    ax = axes;
end


hold(ax,"on")
plot3(ax,nodal_position(1,:),nodal_position(2,:),nodal_position(3,:),"x")
for iNode = 1:num_nodes
    fx = [nodal_position(1,iNode);nodal_force_end(1,iNode)];
    fy = [nodal_position(2,iNode);nodal_force_end(2,iNode)];
    fz = [nodal_position(3,iNode);nodal_force_end(3,iNode)];
    plot3(ax,fx,fy,fz,"-")

    mx = [nodal_position(1,iNode);nodal_moment_end(1,iNode)];
    my = [nodal_position(2,iNode);nodal_moment_end(2,iNode)];
    mz = [nodal_position(3,iNode);nodal_moment_end(3,iNode)];
    plot3(ax,mx,my,mz,"-")
end

daspect(ax,[1,1,1])
xlabel(ax,"x")
ylabel(ax,"y")
zlabel(ax,"z")
end