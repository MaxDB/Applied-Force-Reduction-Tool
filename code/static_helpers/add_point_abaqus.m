function [r,theta,f,E,additional_data] = ...
    add_point_abaqus(applied_force,max_inc,add_data_type,Model,job_id)

JOB_NAME = "static_analysis";
NUM_DIMENSIONS = 6;
project_path = get_project_path;

setup_time_start = tic;

all_dofs = Model.node_mapping(end,1);
num_loadcases = size(applied_force,2);

static_settings = zeros(1,4);
Static_Opts = Model.Static_Options;
static_settings(1) = Static_Opts.initial_time_increment;
static_settings(2) = Static_Opts.total_step_time;
static_settings(3) = Static_Opts.minimum_time_increment;
static_settings(4) = Static_Opts.maximum_time_increment;

new_job = JOB_NAME + "_" + job_id;

%-------------------------------------------------------------------------%
%Open Template
t_id = fopen(project_path + "\fe_templates\abaqus\static_step.inp");
static_template=textscan(t_id,'%s','delimiter','\n');
fclose(t_id);
static_template = static_template{1,1};


for iLine = 1:length(static_template)
    if strfind(static_template{iLine,1},'*Step, name=Static_Step-STEP_NUM, nlgeom=YES, extrapolation=NO, inc=INC_HERE')
        step_def_line = iLine;
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
        static_template{iLine,1} = "*NODE PRINT,SUMMARY=NO,FREQUENCY = " + max_inc;
    end

    if strfind(static_template{iLine,1},'*ENERGY PRINT, FREQUENCY = INC_HERE')
        static_template{iLine,1} = "*ENERGY PRINT, FREQUENCY = " + max_inc;
    end
end


%-------------------------------------------------------------------------%

switch add_data_type
    case "stiffness"
        %Open Template
        t_id = fopen(project_path + "\fe_templates\abaqus\stiffness_output.inp");
        stiffness_template=textscan(t_id,'%s','delimiter','\n');
        fclose(t_id);
        stiffness_template = stiffness_template{1,1};

        for iLine = 1:length(stiffness_template)
            if strfind(stiffness_template{iLine,1},'*Step,name=export_stiffness-STEP_NUM')
                stiffness_def_line = iLine;
            end
        end

    case "perturbation"
        %Open Template
        t_id = fopen(project_path + "\fe_templates\abaqus\static_perturbation_step.inp");
        perturbation_template = textscan(t_id,'%s','delimiter','\n');
        fclose(t_id);
        perturbation_template = perturbation_template{1,1};

        for iLine = 1:length(perturbation_template)
            if strfind(perturbation_template{iLine,1},'*Load Case, Name=')
                loadcase_start_line = iLine;
            end

            if strfind(perturbation_template{iLine,1},'*End Load Case')
                loadcase_end_line = iLine;
            end
        end

        loadcase_template = perturbation_template(loadcase_start_line:loadcase_end_line,:);

        for iLine = 1:length(loadcase_template)
            if strfind(loadcase_template{iLine,1},'LOAD_HERE')
                loadcase_load_line = iLine;
            end
        end


        L_evecs = Model.low_frequency_eigenvectors;
        r_evecs = Model.reduced_eigenvectors;
        h_evecs = [r_evecs,L_evecs];
        num_h_modes = size(h_evecs,2);
        validating_force = (Model.mass*h_evecs)*Static_Opts.perturbation_scale_factor;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
force_transform = Model.mass*Model.reduced_eigenvectors;
coordinate_index = ((1:num_nodes)-1)*NUM_DIMENSIONS;


total_static_steps = num_loadcases;
modal_force = zeros(length(Model.reduced_modes),total_static_steps);

switch add_data_type
    case "none"
        total_steps = total_static_steps;
    case "stiffness"
        total_steps = total_static_steps*2;
    case "perturbation"
        total_steps = total_static_steps*2;

        %perturbation force template
        physical_perturbation_force_bc = validating_force;
        physical_perturbation_force = zeros(all_dofs,num_h_modes);
        physical_perturbation_force(Model.node_mapping(:,1),:) = physical_perturbation_force_bc(Model.node_mapping(:,2),:);

        perturbation_steps = perturbation_template(1:(loadcase_start_line-1));
        loadcase_name_line = zeros(num_h_modes,1);
        for iMode = 1:num_h_modes
            perturbation_force = physical_perturbation_force(:,iMode);

            perturbation_force_label = strings(all_dofs,1);
            for iDimension = 1:NUM_DIMENSIONS
                dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
                perturbation_force_label(dimension_span,1) = force_label(dimension_span,1) + perturbation_force(coordinate_index+iDimension,1);
            end

            loadcase_step = [loadcase_template(1:(loadcase_load_line-1));perturbation_force_label;loadcase_template((loadcase_load_line+1):end)];
            loadcase_name_line(iMode) = length(perturbation_steps) + 1;
            perturbation_steps = [perturbation_steps;loadcase_step];
        end
        perturbation_steps = [perturbation_steps;perturbation_template((loadcase_end_line+1):end)];
end

step_type = strings(total_steps,1);


load_step_counter = 0;
total_step_counter = 0;


try
    input_ID = fopen("temp\" + new_job + ".inp","w");
    fprintf(input_ID,'%s\r\n',geometry{:,1});

    for iLoad = 1:num_loadcases
        load_step_counter = load_step_counter + 1;
        total_step_counter = total_step_counter + 1;

        modal_force(:,load_step_counter) = applied_force(:,iLoad);

        step_force_bc = force_transform*applied_force(:,iLoad);
        step_force = zeros(all_dofs,1);
        step_force(Model.node_mapping(:,1),:) = step_force_bc(Model.node_mapping(:,2),:);
        step_force_label = strings(all_dofs,1);
        for iDimension = 1:NUM_DIMENSIONS
            dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
            step_force_label(dimension_span,1) = force_label(dimension_span,1) + step_force(coordinate_index+iDimension,1);
        end


        static_step = static_template;
        static_step{step_def_line,1} = "*Step, name=STATIC_STEP_" + load_step_counter + ", nlgeom=YES, extrapolation=NO, inc=" + max_inc;

        


        fprintf(input_ID,'%s\r\n',static_step{1:(load_def_line-1),1});
        fprintf(input_ID,'%s\r\n',step_force_label(:));
        fprintf(input_ID,'%s\r\n',static_step{(load_def_line+1):end,1});

        step_type(total_step_counter,1) = "static";

        switch add_data_type
            case "stiffness"
                total_step_counter = total_step_counter + 1;
                stiffness_step = stiffness_template;
                stiffness_step{stiffness_def_line,1} = "*Step,name=STIFFNESS_STEP_" + load_step_counter;

                fprintf(input_ID,'%s\r\n',stiffness_step{:,1});
                step_type(total_step_counter,1) = "stiffness";

            case "perturbation"
                total_step_counter = total_step_counter + 1;
                perturbation_step = perturbation_steps;
                for iMode = 1:num_h_modes
                    name_line = loadcase_name_line(iMode);
                    perturbation_step(name_line) = perturbation_step(name_line) + "step_" +load_step_counter + "-L_" + iMode;
                end

                fprintf(input_ID,'%s\r\n',perturbation_step{:,1});
                step_type(total_step_counter,1) = "perturbation";
        end
    end


catch caught_error
    fclose(input_ID); %ensures file is always closed
    rethrow(caught_error)
end
fclose(input_ID);

setup_time = toc(setup_time_start);
log_message = sprintf("job " + job_id + ": Static input file created: %.1f seconds" ,setup_time);
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
log_message = sprintf("job " + job_id + ": Abaqus static analysis complete: %.1f seconds" ,abaqus_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_processing_time_start = tic;

[displacement,E,additional_data,additional_data_time] = read_abaqus_static_data(new_job,step_type,num_nodes);

displacement_bc = displacement(Model.node_mapping(:,1),:);
disp_transform = force_transform';

r = disp_transform*displacement_bc;
% theta = displacement_bc - Model.reduced_eigenvectors*r;
theta = displacement_bc;

f = modal_force;
step_list = 1:length(E);


data_processing_time = toc(data_processing_time_start) - additional_data_time;
log_message = sprintf("job " + job_id + ": Static data processed: %.1f seconds" ,data_processing_time);
logger(log_message,3)

additional_data_processing_time_start = tic;

switch add_data_type
    case "none"
        additional_data = [];

    case "stiffness"
        additional_data = parse_stiffness(step_list,new_job,Model.num_dof);

    case "perturbation"
        additional_data = additional_data(Model.node_mapping(:,1),:,:);
end
additional_data_processing_time = toc(additional_data_processing_time_start) + additional_data_time;
if add_data_type ~= "none"
    log_message = sprintf("job " + job_id + ": " + add_data_type + " data processed: %.1f seconds" ,additional_data_processing_time);
    logger(log_message,3)
end

%--------------------------

end