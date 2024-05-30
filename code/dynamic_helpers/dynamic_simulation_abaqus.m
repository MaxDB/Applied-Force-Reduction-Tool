function [time,x,x_dot] = dynamic_simulation_abaqus(x_0,x_dot_0,f_r_0,period,min_incs,Model,job_id)
JOB_NAME = "dynamic_analysis";
NUM_DIMENSIONS = 6;
MAX_DYNAMIC_INC = 1e6;
project_path = get_project_path;


setup_time_start = tic;

all_dofs = Model.node_mapping(end,1);
num_dofs = Model.num_dof;

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
        static_template{iLine,1} = "*NODE PRINT,SUMMARY=NO,FREQUENCY = " + max_static_inc;
    end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%Open Template
t_id = fopen(project_path + "\fe_templates\abaqus\dynamic_step.inp");
dynamic_template=textscan(t_id,'%s','delimiter','\n');
fclose(t_id);
dynamic_template = dynamic_template{1,1};


for iLine = 1:length(dynamic_template)
    if strfind(dynamic_template{iLine,1},'*Step, name=Dynamic, nlgeom=YES, inc= INC_HERE')
        dynamic_template{iLine,1} = "*Step, name=Dynamic, nlgeom=YES, inc= " + MAX_DYNAMIC_INC;
    end

    if strfind(dynamic_template{iLine,1},'DYNAMIC_SETTINGS_HERE')
        dynamic_template{iLine,1} = period/1000 + "," + period + "," + period/MAX_DYNAMIC_INC + "," + period/min_incs;
    end
end


%-------------------------------------------------------------------------%
G_ID =fopen("geometry\" + Model.system_name+ "\" + Model.system_name + ".inp");
geometry = textscan(G_ID,'%s','delimiter','\n');
fclose(G_ID);
geometry = geometry{1,1};

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
num_nodes = (all_dofs/NUM_DIMENSIONS);
for iNode = 1:num_nodes
    node_label = instance_name + "." + iNode;
    force_label(iNode) = node_label;

    for iDimension = 1:NUM_DIMENSIONS
        dimension_label = "," + iDimension + ",";
        force_label(iNode+num_nodes*(iDimension-1),1) = node_label + dimension_label;
    end
end

%-------------------------------------------------------------------------%
force_transform = Model.mass*Model.reduced_eigenvectors;
coordinate_index = ((1:num_nodes)-1)*NUM_DIMENSIONS;

total_steps = size(f_r_0,2)*2;

%-------------------------------------------------------------------------%
input_id = fopen("temp\" + new_job + ".inp","w");
fprintf(input_id,'%s\r\n',geometry{:,1});

%-------------------------------------------------------------------------%
step_force_bc = force_transform*f_r_0;
step_force = zeros(all_dofs,1);
step_force(Model.node_mapping(:,1),:) = step_force_bc(Model.node_mapping(:,2),:);
step_force_label = strings(all_dofs,1);
for iDimension = 1:NUM_DIMENSIONS
    dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
    step_force_label(dimension_span,1) = force_label(dimension_span,1) + step_force(coordinate_index+iDimension,1);
end
%-------------------------------------------------------------------------%
step_velocity = zeros(all_dofs,1);
step_velocity(Model.node_mapping(:,1),:) = x_dot_0(Model.node_mapping(:,2),:);
step_velocity_label = strings(all_dofs,1);
for iDimension = 1:NUM_DIMENSIONS
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
dynamic_step = dynamic_template;

fprintf(input_id,'%s\r\n',dynamic_step{1:end});
%-------------------------------------------------------------------------%
fclose(input_id);

setup_time = toc(setup_time_start);
log_message = sprintf("job " + job_id + ": Input file created: %.1f seconds" ,setup_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_cpus = Model.Static_Options.num_fe_cpus;

abaqus_time_start = tic;
project_directory = pwd;
cd temp
[status,cmdout] = system("abaqus job=" + new_job + " cpus=" + num_cpus); %#ok<ASGLU>

while ~isfile(new_job + ".dat")
    pause(0.1)
end

while isfile(new_job + ".lck")
    pause(0.1)
end
cd(project_directory)

abaqus_time = toc(abaqus_time_start);
log_message = sprintf("job " + job_id + ": Abaqus dynamic analysis complete: %.1f seconds" ,abaqus_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_start_pattern = "STEP" + whitespacePattern + digitsPattern + whitespacePattern + "INCREMENT" + whitespacePattern + "1";
increment_start_pattern = "INCREMENT" + whitespacePattern + digitsPattern + whitespacePattern + "SUMMARY";
disp_table_pattern = "NODE FOOT-" + whitespacePattern +  "U1";
vel_table_pattern = "NODE FOOT-" + whitespacePattern +  "V1";
increment_time_pattern = "TIME INCREMENT COMPLETED";



%------------------------------------------------------------------------%
dat_ID = fopen("temp\" + new_job + ".dat");
abaqus_data = textscan(dat_ID,'%s','delimiter','\n');
fclose(dat_ID);
abaqus_data = abaqus_data{1,1};
%------------------------------------------------------------------------%
step_start_lines = [find(matches(abaqus_data,step_start_pattern,'IgnoreCase',true),2);length(abaqus_data)];
inc_start_lines = [find(matches(abaqus_data,increment_start_pattern,'IgnoreCase',true));length(abaqus_data)];

%Initial displacement
step_span = step_start_lines(1):(step_start_lines(2)-1);
step_data = abaqus_data(step_span,1);
disp_table_start = find(startsWith(step_data,disp_table_pattern,'IgnoreCase',true),1);
disp_table_span = disp_table_start:size(step_span,2);
disp_table_data = step_data(disp_table_span,1);
disp_0 = read_abaqus_table(disp_table_data,num_nodes,NUM_DIMENSIONS);
disp_0_bc = disp_0(Model.node_mapping(:,1),:);


%Dynamic increments
num_increments = size(inc_start_lines,1) - 2;
time = zeros(1,num_increments+1);
displacement = zeros(num_dofs,num_increments);
velocity = zeros(num_dofs,num_increments);

for iInc = 1:num_increments
    inc_span = inc_start_lines(iInc+1):(inc_start_lines(iInc+2)-1);
    inc_data = abaqus_data(inc_span,1);
    increment_time_line = find(startsWith(inc_data,increment_time_pattern,'IgnoreCase',true),1);
    disp_table_start = find(startsWith(inc_data,disp_table_pattern,'IgnoreCase',true),1);
    vel_table_start = find(startsWith(inc_data,vel_table_pattern,'IgnoreCase',true),1);
     
    increment_time_line_data = textscan(inc_data{increment_time_line},"%s %s %s %f");
    time(iInc+1) = increment_time_line_data{1,4} + time(iInc);

    disp_table_span = disp_table_start:(vel_table_start-3);
    disp_table_data = inc_data(disp_table_span,1);
    disp_pre_bc = read_abaqus_table(disp_table_data,num_nodes,NUM_DIMENSIONS);
    displacement(:,iInc) = disp_pre_bc(Model.node_mapping(:,1),:);

    vel_table_span = vel_table_start:size(inc_span,2);
    vel_table_data = inc_data(vel_table_span,1);
    vel_pre_bc = read_abaqus_table(vel_table_data,num_nodes,NUM_DIMENSIONS);
    velocity(:,iInc) = vel_pre_bc(Model.node_mapping(:,1),:);

end


x = [disp_0_bc,displacement];
x_dot = [x_dot_0,velocity];
end



function step_displacement = read_abaqus_table(table_data,num_nodes,num_dimensions)
% node_output_pattern = "N O D E   O U T P U T";
% output_start_line = find(startsWith(table_data,node_output_pattern,'IgnoreCase',true),1);
% node_data = table_data(output_start_line:end,1);

step_displacement = zeros(num_nodes,num_dimensions);
for iLine = 1:length(table_data)
    line = table_data{iLine,1};
    if isempty(line)
        continue
    end

    line_data = textscan(line,"%u %f %f %f %f %f %f");
    if isempty(line_data{1,7})
        continue
    end

    node_num = line_data{1,1};
    step_displacement(node_num,:) = [line_data{1,2:end}];
end

num_dofs = num_nodes*num_dimensions;
step_displacement = reshape(step_displacement',num_dofs,1);
end