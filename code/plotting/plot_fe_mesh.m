function plot_fe_mesh(Model,displacement,colour_map)
%input_file = nargin == 1;?
input_file = 0;

if input_file
    geometry = Model;
else
    geometry = load_geometry(Model);
end
mesh_data = get_mesh_data(geometry);
mesh_data = mesh_data{1};
is_3d = mesh_data.type == "continuous 3D";

node_position = read_abaqus_node_position(geometry);
element_members = read_abaqus_element_membership(geometry);


plot_colour_map = nargin > 1;

if plot_colour_map
    num_dimensions = get_num_node_dimensions(Model);
    num_nodes = size(node_position,1);
    num_dof = num_nodes*num_dimensions;
    if size(displacement,1) > num_nodes
        if size(displacement,1) < num_dof
            displacement_bc = displacement;
            displacement = zeros(num_dof,1);
            node_map = Model.node_mapping;
            displacement(node_map(:,1)) = displacement_bc(node_map(:,2));
        end

        displacement = reshape(displacement,[num_dimensions,num_nodes])';

    end

    node_position = node_position + displacement(:,1:num_dimensions/2);
    if nargin == 2
        colour_map = node_position;
    end
    patch_options = {"FaceVertexCData",colour_map,"FaceColor","interp"};
else
    patch_options = {"FaceColor","green"};

end

figure;
ax = gca;
box on

if ~is_3d
    p = patch(ax,"Faces",element_members,"Vertices",node_position,patch_options{:});
    colorbar
    dt = datatip(p);
    data_tip_row = dataTipTextRow("Node",element_members');
    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
    if plot_colour_map
        data_tip_row = dataTipTextRow("Value",colour_map(element_members'));
        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
    end
    delete(dt);
else
    p = plot3(node_position(:,1),node_position(:,2),node_position(:,3),".");
end






ax.DataAspectRatio = [1,1,1];
end