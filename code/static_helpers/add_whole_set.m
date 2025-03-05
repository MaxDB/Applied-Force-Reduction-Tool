function geometry = add_whole_set(geometry,instance_name)
SINGLE_ASTERIX_START = "^\*(?!\*)";
ELEMENT_LINE  = "*Element";
NODE_LINE = "*Node";
ASSEMBLY_END = "*End Assembly";

SET_NAME = "set-all";

command_start = find(~cellfun(@isempty,(regexp(geometry,SINGLE_ASTERIX_START,"once"))));
commands = geometry(command_start);
% command_start(end) = size(geometry,1) + 1;

node_def_line = find(startsWith(commands,NODE_LINE));
element_def_line = find(startsWith(commands,ELEMENT_LINE));

node_def_ends = command_start(node_def_line + 1) - 1;
element_def_ends = command_start(element_def_line + 1) - 1;

max_node_line = geometry{node_def_ends};
max_node_line_parts = split(max_node_line,",");
max_node = max_node_line_parts{1};

max_element_line = geometry{element_def_ends};
max_element_line_parts = split(max_element_line,",");
max_element = max_element_line_parts{1};

nset_def = "*Nset, nset=" + SET_NAME + ", instance=" + instance_name + ", generate";
nset_span = "1, " + max_node + ", 1";

elset_def = "*Elset, elset=" + SET_NAME + ", instance=" + instance_name + ", generate";
elset_span = "1, " + max_element + ", 1";

all_set_def = {nset_def;nset_span;elset_def;elset_span};
all_set_def = cellfun(@convertStringsToChars,all_set_def,"UniformOutput",false);

assembly_end_line = command_start((startsWith(commands,ASSEMBLY_END)));
geometry = [geometry(1:(assembly_end_line-1));all_set_def;geometry(assembly_end_line:end)];
end