function [stress,stress_labels,section_points] = stress_simulation_abaqus(x_0,x_dot_0,f_r_0,Model,job_id)
JOB_NAME = "stress_analysis";
num_dimensions = get_num_node_dimensions(Model);
project_path = get_project_path;


setup_time_start = tic;

all_dofs = Model.node_mapping(end,1);

static_settings = zeros(1,4);
Static_Opts = Model.Static_Options;
static_settings(1) = Static_Opts.initial_time_increment;
static_settings(2) = Static_Opts.total_step_time;
static_settings(3) = Static_Opts.minimum_time_increment;
static_settings(4) = Static_Opts.maximum_time_increment;
max_static_inc = Static_Opts.maximum_step_increments*Static_Opts.num_loadcases;

new_job = JOB_NAME + "_" + job_id;



%-------------------------------------------------------------------------%
%Open Template
t_id = fopen(project_path + "\fe_templates\abaqus\static_ic_step.inp");
static_template=textscan(t_id,'%s','delimiter','\n');
fclose(t_id);
static_template = static_template{1,1};


for iLine = 1:length(static_template)
    if strfind(static_template{iLine,1},'VELOCITY_HERE')
        velocity_def_line = iLine;
    end

    if strfind(static_template{iLine,1},'*Step, name=Static_Step-STEP_NUM, nlgeom=YES, extrapolation=NO, inc=INC_HERE')
        static_template{iLine,1} = "*Step, name=Static_Step, nlgeom=YES, extrapolation=NO, inc=" + max_static_inc;
    end

    if strfind(static_template{iLine,1},'SETTINGS_HERE')
        settings_string = "";
        for iSetting = 1:length(static_settings)
            settings_string = settings_string + static_settings(iSetting) + ","; %#ok<*AGROW>
        end
        static_template{iLine,1} = settings_string;
    end

    if strfind(static_template{iLine,1},'LOAD_HERE')
        load_def_line = iLine;
        static_step{load_def_line-1} = "*Cload, OP = NEW";
    end

    if strfind(static_template{iLine,1},'*NODE PRINT,SUMMARY=NO,FREQUENCY = INC_HERE')
        static_template{iLine,1} = "*EL PRINT,SUMMARY=NO, FREQUENCY = " +  max_static_inc + ", POSITION=AVERAGED AT NODES";
        static_template{iLine+1,1} = "S";
    end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
geometry = load_geometry(Model);

for iLine = 1:length(geometry)
    if strfind(geometry{iLine,1},"*Instance, name=")
        instance_def = erase(geometry{iLine,1},"*Instance, name=");
        instance_def = split(instance_def,",");
        instance_name = instance_def{1,1};
    end
end


%-------------------------------------------------------------------------%
%%% Create forcing tempate
force_label = strings(all_dofs,1);
num_nodes = (all_dofs/num_dimensions);
for iNode = 1:num_nodes
    node_label = instance_name + "." + iNode;
    force_label(iNode) = node_label;

    for iDimension = 1:num_dimensions
        dimension_label = "," + iDimension + ",";
        force_label(iNode+num_nodes*(iDimension-1),1) = node_label + dimension_label;
    end
end

%-------------------------------------------------------------------------%
force_transform = Model.mass*Model.reduced_eigenvectors;
coordinate_index = ((1:num_nodes)-1)*num_dimensions;

%-------------------------------------------------------------------------%
input_id = fopen("temp\" + new_job + ".inp","w");
fprintf(input_id,'%s\r\n',geometry{:,1});

%-------------------------------------------------------------------------%
step_force_bc = force_transform*f_r_0;
step_force = zeros(all_dofs,1);
step_force(Model.node_mapping(:,1),:) = step_force_bc(Model.node_mapping(:,2),:);
step_force_label = strings(all_dofs,1);
for iDimension = 1:num_dimensions
    dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
    step_force_label(dimension_span,1) = force_label(dimension_span,1) + step_force(coordinate_index+iDimension,1);
end
%-------------------------------------------------------------------------%
step_velocity = zeros(all_dofs,1);
step_velocity(Model.node_mapping(:,1),:) = x_dot_0(Model.node_mapping(:,2),:);
step_velocity_label = strings(all_dofs,1);
for iDimension = 1:num_dimensions
    dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
    step_velocity_label(dimension_span,1) = force_label(dimension_span,1) + step_velocity(coordinate_index+iDimension,1);
end

%-------------------------------------------------------------------------%
static_step = static_template;

fprintf(input_id,'%s\r\n',static_step{1:(velocity_def_line-1),1});
fprintf(input_id,'%s\r\n',step_velocity_label(:));
fprintf(input_id,'%s\r\n',static_step{(velocity_def_line+1):(load_def_line-1),1});
fprintf(input_id,'%s\r\n',step_force_label(:));
fprintf(input_id,'%s\r\n',static_step{(load_def_line+1):end,1});


%-------------------------------------------------------------------------%
fclose(input_id);

setup_time = toc(setup_time_start);
log_message = sprintf("job " + job_id + ": Input file created: %.1f seconds" ,setup_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_cpus = Model.Static_Options.num_fe_cpus;

abaqus_time_start = tic;
[status,cmdout] = run_abaqus_job(new_job,"num_cpus",num_cpus,"interactive",3); %#ok<ASGLU>

abaqus_time = toc(abaqus_time_start);
log_message = sprintf("job " + job_id + ": Abaqus dynamic analysis complete: %.1f seconds" ,abaqus_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
element_table_pattern = "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS WITH TYPE " + alphanumericsPattern + " AVERAGED AT THE NODES";



%------------------------------------------------------------------------%
dat_ID = fopen("temp\" + new_job + ".dat");
abaqus_data = textscan(dat_ID,'%s','delimiter','\n');
fclose(dat_ID);
abaqus_data = abaqus_data{1,1};
%------------------------------------------------------------------------%
table_lines = find(matches(abaqus_data,element_table_pattern,'IgnoreCase',true));
num_tables = size(table_lines,1);

for iTable = 1:num_tables
    table_def_line = abaqus_data{table_lines(iTable),1};
    if contains(table_def_line,"SPRING","IgnoreCase",true)
        continue
    end

    table_range = [0;0];
    table_lines_max = (length(abaqus_data)-table_lines(iTable));
    for iLine = 1:table_lines_max
        line = abaqus_data{table_lines(iTable)+iLine,1};
        if isempty(line) || startsWith(line,lettersPattern)
            if table_range(1) && ~table_range(2)
                table_range(2) = table_lines(iTable)+(iLine-1);
            end
            continue
        end
        if ~table_range(1)
            table_range(1) = table_lines(iTable)+iLine;
        end
    end

    stress_table = abaqus_data(table_range(1):table_range(2));
    stress_cell = cellfun(@(line) textscan(line,"%f"),stress_table);
    first_node = stress_cell{1,1}(1);
    section_point_lines = startsWith(stress_table,string(first_node) + whitespacePattern);
    section_points = cellfun(@(line) line(2), stress_cell(section_point_lines))';
    num_section_points = size(section_points,1);

    output_line = split(abaqus_data(table_range(1)-3));
    stress_label_indices = matches(output_line,"S" + digitsPattern);
    stress_labels = output_line(stress_label_indices);
    stress_labels = cellfun(@(line) string(line),stress_labels)';
    num_stress_outputs = size(stress_labels,2);

    num_rows = size(stress_table,1);
    stress = zeros(num_nodes,num_stress_outputs,num_section_points);
    for iRow = 1:num_rows
        line = stress_cell{iRow,1};
        node = line(1);
        section_point = line(2);
        principal_stress = line(3:end);
        stress(node,:,section_point == section_points) = principal_stress;
    end
    

end

end