function element_members = read_abaqus_element_membership(geometry)

SECTION_START_PATTERN = "*";
ELEMENT_START_PATTERN = "*Element";
section_starts = find(startsWith(geometry,SECTION_START_PATTERN,'IgnoreCase',true));

%----- Elements
element_start_index = find(startsWith(geometry(section_starts),ELEMENT_START_PATTERN,'IgnoreCase',true),1);
element_range = section_starts(element_start_index + [0,1]);
num_elements = diff(element_range) - 1;


%check element def line height
for iElement = 1:num_elements
    line = geometry{element_range(1) + iElement};
    line_data = textscan(line,"%u", "Delimiter",",");
    line_data = line_data{1};
    if line_data(1) == 2
        element_def_height = iElement - 1;
        break
    end
end

num_elements = num_elements/element_def_height;
line_length = 0;
for iLine = 1:element_def_height
    element_line = geometry{element_range(1)+iLine};
    line_length = line_length + length(split(element_line,","));
end

element_members = zeros(num_elements,line_length-1);

for iElement = 1:num_elements
    element_line_data = zeros(line_length-1,1);
    element_counter = 0;
    for iLine = 1:element_def_height
        line = geometry{element_range(1) + 2*(iElement-1) + iLine};
        line_data = textscan(line,"%u", "Delimiter",",");
        line_data = line_data{1};
        if iLine == 1
            element_id = line_data(1);
            line_data(1) = [];
        end
        element_line_length = size(line_data,1);
        element_span = (element_counter + 1):(element_counter + element_line_length);
        element_line_data(element_span) = line_data;
        element_counter = element_line_length;
    end
    element_members(element_id,:) = element_line_data;
    % element_members(line_data(1),:) = line_data;
end

end