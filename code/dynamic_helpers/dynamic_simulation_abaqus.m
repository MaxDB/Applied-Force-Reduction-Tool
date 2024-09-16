function [time,x,x_dot,energy] = dynamic_simulation_abaqus(x_0,x_dot_0,f_r_0,period,num_periods,min_incs,initial_time,FE_Force_Data,Model,job_id)
JOB_NAME = "dynamic_analysis";
NUM_DIMENSIONS = 6;
MAX_DYNAMIC_INC = 1e6;

AMPLITUDE_TYPE = "periodic_amplitude";
TEMPLATE_PATH = "\fe_templates\abaqus\";
project_path = get_project_path;


setup_time_start = tic;

all_dofs = Model.node_mapping(end,1);
num_dofs = Model.num_dof;
num_nodes = (all_dofs/NUM_DIMENSIONS);

static_settings = zeros(1,4);
Static_Opts = Model.Static_Options;
static_settings(1) = Static_Opts.initial_time_increment;
static_settings(2) = Static_Opts.total_step_time;
static_settings(3) = Static_Opts.minimum_time_increment;
static_settings(4) = Static_Opts.maximum_time_increment;
max_static_inc = Static_Opts.maximum_step_increments*Static_Opts.num_loadcases;


new_job = JOB_NAME + "_" + job_id(1);

restart_write = size(job_id,2) == 2;
restart_read = 0;

if restart_write
    new_job = new_job + "_" + job_id(2);
    restart_step = job_id(2);
    restart_read = restart_step > 1;
end

if restart_read
    old_job = JOB_NAME + "_" + job_id(1) + "_" + (job_id(2)-1);
end

%-------------------------------------------------------------------------%
%Open Template
if ~isempty(FE_Force_Data)
    frequency = 2*pi/period;
    harmonic_coefficients = shift_harmonics(FE_Force_Data.harmonic_coefficients,initial_time,frequency);
    initial_time = 0;
    % harmonic_coefficients = FE_Force_Data.harmonic_coefficients;

    t_id = fopen(project_path + TEMPLATE_PATH + AMPLITUDE_TYPE + ".inp");
    amp_template=textscan(t_id,'%s','delimiter','\n');
    fclose(t_id);
    amp_template = amp_template{1,1};

    for iLine = 1:length(amp_template)
        if strfind(amp_template{iLine,1},'1, freq, t0, A0')
            amp_template{iLine,1} = "1, " + frequency + ", " + initial_time + ", " + harmonic_coefficients(1);
            amp_template{iLine+1,1} = harmonic_coefficients(2) + ", " + harmonic_coefficients(3);
        end
    end
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%Open Template
t_id = fopen(project_path + TEMPLATE_PATH + "dynamic_step.inp");
dynamic_template=textscan(t_id,'%s','delimiter','\n');
fclose(t_id);
dynamic_template = dynamic_template{1,1};


for iLine = 1:length(dynamic_template)
    if strfind(dynamic_template{iLine,1},'*Step, name=Dynamic, nlgeom=YES, inc= INC_HERE')
        dynamic_template{iLine,1} = "*Step, name=Dynamic, nlgeom=YES, inc= " + MAX_DYNAMIC_INC;
    end

    if strfind(dynamic_template{iLine,1},'DYNAMIC_SETTINGS_HERE')
        dynamic_template{iLine,1} = period/1000 + "," + period*num_periods + "," + period/MAX_DYNAMIC_INC + "," + period/min_incs;
    end

    if ~isempty(FE_Force_Data)
        if strfind(dynamic_template{iLine,1},'*CLOAD, OP=NEW')
            dynamic_template{iLine,1} = "*CLOAD, amplitude=" + AMPLITUDE_TYPE +  ", OP=NEW";
        end
        if strfind(dynamic_template{iLine,1},'**DYNAMIC_LOAD_HERE')
            dynamic_load_def_line = iLine;
        end
    end

    if restart_write
        if strfind(dynamic_template{iLine,1},'**Restart, write, frequency = INC_HERE')
            dynamic_template{iLine,1} = "*Restart, write, frequency =" + MAX_DYNAMIC_INC;
        end
    end

    if restart_read
        if strfind(dynamic_template{iLine,1},'**Restart, read, step = STEP_HERE, inc = INC_HERE')
            dynamic_template{iLine,1} = "*Restart, read, step = " + restart_step;
        end
    end
end


%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%Open Template

t_id = fopen(project_path + TEMPLATE_PATH + "static_ic_step.inp");
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

geometry = load_geometry(Model);

for iLine = 1:length(geometry)
    if strfind(geometry{iLine,1},"*Instance, name=")
        instance_def = erase(geometry{iLine,1},"*Instance, name=");
        instance_def = split(instance_def,",");
        instance_name = instance_def{1,1};
    end


    if strfind(geometry{iLine,1},"*Density")
        density_def_line = iLine;
    end

end
if ~isempty(FE_Force_Data)
    alpha = FE_Force_Data.alpha;
    beta = FE_Force_Data.beta;
    damping_def = "*Damping, alpha = " + alpha + ", beta = " + beta;
    geometry = [geometry(1:(density_def_line-1));damping_def;geometry(density_def_line:end);amp_template];
end

%-------------------------------------------------------------------------%
%%% Create forcing tempate
force_label = strings(all_dofs,1);
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


%-------------------------------------------------------------------------%
%create job
input_id = fopen("temp\" + new_job + ".inp","w");


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
if ~isempty(FE_Force_Data)
    % dyn_force_bc = FE_Force_Data.amplitude*FE_Force_Data.force_shape;
    dyn_force_bc = FE_Force_Data.amplitude*Model.mass*FE_Force_Data.force_shape;
    dyn_force = zeros(all_dofs,1);
    dyn_force(Model.node_mapping(:,1),:) = dyn_force_bc(Model.node_mapping(:,2),:);
    dyn_force_label = strings(all_dofs,1);
    for iDimension = 1:NUM_DIMENSIONS
        dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
        dyn_force_label(dimension_span,1) = force_label(dimension_span,1) + dyn_force(coordinate_index+iDimension,1);
    end
end
%-------------------------------------------------------------------------%
if ~restart_read
    fprintf(input_id,'%s\r\n',geometry{:,1});
    static_step = static_template;

    fprintf(input_id,'%s\r\n',static_step{1:(velocity_def_line-1),1});
    fprintf(input_id,'%s\r\n',step_velocity_label(:));
    fprintf(input_id,'%s\r\n',static_step{(velocity_def_line+1):(load_def_line-1),1});
    fprintf(input_id,'%s\r\n',step_force_label(:));
    fprintf(input_id,'%s\r\n',static_step{(load_def_line+1):end,1});
end
%-------------------------------------------------------------------------%
dynamic_step = dynamic_template;
if ~isempty(FE_Force_Data)
    fprintf(input_id,'%s\r\n',dynamic_step{1:(dynamic_load_def_line-1)});
    fprintf(input_id,'%s\r\n',dyn_force_label(:));
    fprintf(input_id,'%s\r\n',dynamic_step{(dynamic_load_def_line+1):end,1});
else
    fprintf(input_id,'%s\r\n',dynamic_step{:});
end



%-------------------------------------------------------------------------%
fclose(input_id);

setup_time = toc(setup_time_start);
log_message = sprintf("job " + job_id(1) + ": Input file created: %.1f seconds" ,setup_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_cpus = Model.Static_Options.num_fe_cpus;

abaqus_time_start = tic;
project_directory = pwd;
cd temp
if ~restart_read
    [status,cmdout] = system("abaqus job=" + new_job + " cpus=" + num_cpus); %#ok<ASGLU>
else
    [status,cmdout] = system("abaqus job=" + new_job + " oldjob=" + old_job + " cpus=" + num_cpus); %#ok<ASGLU>
end

while ~isfile(new_job + ".dat")
    pause(0.1)
end

while isfile(new_job + ".lck")
    pause(0.1)
end
cd(project_directory)

abaqus_time = toc(abaqus_time_start);
log_message = sprintf("job " + job_id(1) + ": Abaqus dynamic analysis complete: %.1f seconds" ,abaqus_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_processing_time_start = tic;

step_start_pattern = "STEP" + whitespacePattern + digitsPattern + whitespacePattern + "INCREMENT" + whitespacePattern + "1";
increment_start_pattern = "INCREMENT" + whitespacePattern + digitsPattern + whitespacePattern + "SUMMARY";
disp_table_pattern = "NODE FOOT-" + whitespacePattern +  "U1";
vel_table_pattern = "NODE FOOT-" + whitespacePattern +  "V1";
increment_time_pattern = "TIME INCREMENT COMPLETED";

potential_pattern = "RECOVERABLE STRAIN ENERGY";
kinetic_pattern = "KINETIC ENERGY";
work_pattern = "EXTERNAL WORK";
dissipation_pattern = "VISCOUS DISSIPATION";

%------------------------------------------------------------------------%
dat_ID = fopen("temp\" + new_job + ".dat");
abaqus_data = textscan(dat_ID,'%s','delimiter','\n');
fclose(dat_ID);
abaqus_data = abaqus_data{1,1};
%------------------------------------------------------------------------%
step_start_lines = [find(matches(abaqus_data,step_start_pattern,'IgnoreCase',true),2);length(abaqus_data)];
inc_start_lines = [find(matches(abaqus_data,increment_start_pattern,'IgnoreCase',true));length(abaqus_data)];

%Initial displacement
if ~restart_read
    step_span = step_start_lines(1):(step_start_lines(2)-1);
    step_data = abaqus_data(step_span,1);
    disp_table_start = find(startsWith(step_data,disp_table_pattern,'IgnoreCase',true),1);
    disp_table_span = disp_table_start:size(step_span,2);
    disp_table_data = step_data(disp_table_span,1);
    disp_0 = read_abaqus_table(disp_table_data,num_nodes,NUM_DIMENSIONS);
    disp_0_bc = disp_0(Model.node_mapping(:,1),:);
end

%Dynamic increments
num_increments = size(inc_start_lines,1) - 2;
time = zeros(1,num_increments+1);
displacement = zeros(num_dofs,num_increments);
velocity = zeros(num_dofs,num_increments);
potential_energy = zeros(1,num_increments);
kinetic_energy = zeros(1,num_increments);
external_work = zeros(1,num_increments);
dissipated_energy = zeros(1,num_increments);

for iInc = 1:num_increments
    inc_span = inc_start_lines(iInc+1):(inc_start_lines(iInc+2)-1);
    inc_data = abaqus_data(inc_span,1);
    increment_time_line = find(startsWith(inc_data,increment_time_pattern,'IgnoreCase',true),1);
    disp_table_start = find(startsWith(inc_data,disp_table_pattern,'IgnoreCase',true),1);
    vel_table_start = find(startsWith(inc_data,vel_table_pattern,'IgnoreCase',true),1);

    increment_time_line_data = textscan(inc_data{increment_time_line},"%s %s %s %f");
    time(iInc+1) = increment_time_line_data{1,4} + time(iInc);

    potential_line_def = find(startsWith(inc_data,potential_pattern,'IgnoreCase',true),1);
    potential_line = textscan(inc_data{potential_line_def},"%s %s %s %f");
    potential_energy(:,iInc) = potential_line{1,end};

    kinetic_line_def = find(startsWith(inc_data,kinetic_pattern,'IgnoreCase',true),1);
    kinetic_line = textscan(inc_data{kinetic_line_def},"%s %s %f");
    kinetic_energy(:,iInc) = kinetic_line{1,end};

    external_work_line_def = find(startsWith(inc_data,work_pattern,'IgnoreCase',true),1);
    external_work_line = textscan(inc_data{external_work_line_def},"%s %s %f");
    external_work(:,iInc) = external_work_line{1,end};

    dissipated_energy_line_def = find(startsWith(inc_data,dissipation_pattern,'IgnoreCase',true),1);
    dissipated_energy_line = textscan(inc_data{dissipated_energy_line_def},"%s %s %s %s %s %f");
    dissipated_energy(:,iInc) = dissipated_energy_line{1,end};


    disp_table_span = disp_table_start:(vel_table_start-3);
    disp_table_data = inc_data(disp_table_span,1);
    disp_pre_bc = read_abaqus_table(disp_table_data,num_nodes,NUM_DIMENSIONS);
    displacement(:,iInc) = disp_pre_bc(Model.node_mapping(:,1),:);

    vel_table_span = vel_table_start:size(inc_span,2);
    vel_table_data = inc_data(vel_table_span,1);
    vel_pre_bc = read_abaqus_table(vel_table_data,num_nodes,NUM_DIMENSIONS);
    velocity(:,iInc) = vel_pre_bc(Model.node_mapping(:,1),:);

end

energy.potential = potential_energy;
energy.work = external_work;
energy.kinetic = kinetic_energy;
energy.dissipated = dissipated_energy;

if ~restart_read
    x = [disp_0_bc,displacement];
    x_dot = [x_dot_0,velocity];
else
    x = displacement;
    x_dot = velocity;
    time = time(2:end);
end

data_processing_time = toc(data_processing_time_start);
log_message = sprintf("job " + job_id(1) + ": Dynamic data processed: %.1f seconds" ,data_processing_time);
logger(log_message,3)
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