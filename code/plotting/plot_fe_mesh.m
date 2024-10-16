function plot_fe_mesh(Model,displacement,colour_map)
NUM_DIMENSIONS = 6;

geometry = load_geometry(Model);
node_position = read_abaqus_node_position(geometry);
element_members = read_abaqus_element_membership(geometry);

num_nodes = size(node_position,1);
num_dof = num_nodes*NUM_DIMENSIONS;
if size(displacement,1) > num_nodes
    if size(displacement,1) < num_dof
        displacement_bc = displacement;
        displacement = zeros(num_dof,1);
        node_map = Model.node_mapping;
        displacement(node_map(:,1)) = displacement_bc(node_map(:,2));
    end

    displacement = reshape(displacement,[NUM_DIMENSIONS,num_nodes])';

end

plot_colour_map = nargin > 1;


if plot_colour_map
    node_position = node_position + displacement(:,1:NUM_DIMENSIONS/2);
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

ax.DataAspectRatio = [1,1,1];
end