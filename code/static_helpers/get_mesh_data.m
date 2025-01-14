function mesh_data = get_mesh_data(geometry)
ELEMENT_DEF = "*Element";
ASSEMBLY_DEF = "*Assembly";
TYPE_DEF = 'type=';

assembly_def_line = find(startsWith(geometry,ASSEMBLY_DEF,'IgnoreCase',true),1,"last");
element_def_lines = find(startsWith(geometry(1:assembly_def_line),ELEMENT_DEF,'IgnoreCase',true));
elements_def = geometry(element_def_lines);

num_element_types = length(element_def_lines);
if num_element_types > 1
    warning("multiple element types detected")
end

mesh_data = cell(num_element_types,1);
for iElement = 1:num_element_types

    element_def_line = elements_def{iElement};
    element_def_parts = split(element_def_line,",");
    element_type_part = element_def_parts{2};
    type_start = strfind(element_type_part,TYPE_DEF);
    element_def = element_type_part((type_start+length(TYPE_DEF)):end);



    switch element_def(1)
        case 'C'
            element_type = "continuous 3D";
            element_dimension = 3;
        case 'B'
            element_type = "beam";
            element_dimension = 6;
        case 'S'
            element_type = "shell";
            element_dimension = 6;
        otherwise
            error("Unsupported element type: " + element_def)
    end

end

clear("Element_Data")
Element_Data.definition = element_def;
Element_Data.type = element_type;
Element_Data.dimension = element_dimension;
mesh_data{iElement} = Element_Data;

end