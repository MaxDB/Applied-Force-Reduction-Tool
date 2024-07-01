function element_members = read_abaqus_element_membership(geometry)

SECTION_START_PATTERN = "*";
ELEMENT_START_PATTERN = "*Element";
section_starts = find(startsWith(geometry,SECTION_START_PATTERN,'IgnoreCase',true));

%----- Elements
element_start_index = find(startsWith(geometry(section_starts),ELEMENT_START_PATTERN,'IgnoreCase',true),1);
element_range = section_starts(element_start_index + [0,1]);
num_elements = diff(element_range) - 1;

element_line = geometry{element_range(1)+1};
line_length = length(split(element_line,","));

element_members = zeros(num_elements,line_length-1);
for iElement = 1:num_elements
    line = geometry{element_range(1) + iElement};
    line_data = textscan(line,"%u", "Delimiter",",");
    line_data = line_data{1};
    element_members(line_data(1),:) = line_data(2:end);
end

end