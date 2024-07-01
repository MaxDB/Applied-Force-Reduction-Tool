function node_position = read_abaqus_node_position(geometry)
NUM_LINEAR_DIMENSIONS = 3;

SECTION_START_PATTERN = "*";
NODE_START_PATTERN = "*Node";
section_starts = find(startsWith(geometry,SECTION_START_PATTERN,'IgnoreCase',true));

%----- Nodes
node_start_index = find(startsWith(geometry(section_starts),NODE_START_PATTERN,'IgnoreCase',true),1);
node_range = section_starts(node_start_index + [0,1]);

num_nodes = diff(node_range) - 1;

node_position = zeros(num_nodes,NUM_LINEAR_DIMENSIONS);
for iNode = 1:num_nodes
    line = geometry{node_range(1) + iNode};
    line_data = textscan(line,"%u %f %f %f","Delimiter",",");
    node_position(line_data{1},:) = [line_data{2:end}];
end

end