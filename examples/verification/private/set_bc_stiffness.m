function set_bc_stiffness(system_name,x_spring_stiffness,y_spring_stiffness)
STIFFNESS_PATTERN = "STIFFNESS_HERE";

geometry_path = "geometry\" + system_name + "\";
geometry_file_path = geometry_path + system_name + ".inp";
geometry_template_path = geometry_path + system_name + "_template.inp";

gt_id =fopen(geometry_template_path);
geometry_template = textscan(gt_id,'%s','delimiter','\n');
fclose(gt_id);
geometry_template = geometry_template{1,1};

stiffness_lines = find(endsWith(geometry_template,STIFFNESS_PATTERN));
x_line = stiffness_lines(startsWith(geometry_template(stiffness_lines),"X"));
y_line = stiffness_lines(startsWith(geometry_template(stiffness_lines),"Y"));

if isfile(geometry_file_path)
    g_id =fopen(geometry_file_path);
    geometry = textscan(g_id,'%s','delimiter','\n');
    fclose(g_id);
    geometry = geometry{1,1};
    if geometry{x_line} == string(x_spring_stiffness) && geometry{y_line} == string(y_spring_stiffness)
        return
    end
end

geometry_template{x_line} = string(x_spring_stiffness);
geometry_template{y_line} = string(y_spring_stiffness);


reset_geometry_directory(system_name)

try
    input_ID = fopen(geometry_file_path,"w");
    fprintf(input_ID,'%s\r\n',geometry_template{:,1});
catch caught_error
    fclose(input_ID); %ensures file is always closed
    rethrow(caught_error)
end
fclose(input_ID);
end