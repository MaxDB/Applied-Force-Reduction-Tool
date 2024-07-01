function plot_fe_mesh(Model,displacement,colour_map)
NUM_DIMENSIONS = 6;

geometry = load_geometry(Model);
node_position = read_abaqus_node_position(geometry);
element_members = read_abaqus_element_membership(geometry);


plot_colour_map = nargin > 1;

if plot_colour_map
    node_position = node_position + displacement(:,1:NUM_DIMENSIONS/2);
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