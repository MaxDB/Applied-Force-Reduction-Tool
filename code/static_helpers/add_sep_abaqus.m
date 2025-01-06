function [r,theta,f,E,additional_data,sep_id] = ...
    add_sep_abaqus(force_ratio,num_loadcases,Static_Opts,max_inc,add_data_type,clean_data,Model,job_id,initial_load,restart_type)

JOB_NAME = "static_analysis";
num_dimensions = get_num_node_dimensions(Model);
project_path = get_project_path;

setup_time_start = tic;

all_dofs = Model.node_mapping(end,1);
% deactivated_dofs = mod(all_dofs,NUM_DIMENSIONS) ~= 0
deactivated_dofs = 0;
if deactivated_dofs
    % replace with real solution
    EXCLUDED_DOF = 242*[4,5,6]; %#ok<UNRCH>
    all_dofs = floor(all_dofs/num_dimensions)*num_dimensions +num_dimensions;
end
num_seps = size(force_ratio,2);

if ~exist("restart_type","var")
    restart_type = zeros(1,num_seps);
end

if isscalar(num_loadcases)
    num_loadcases = ones(1,num_seps)*num_loadcases;
end

if ~exist("initial_load","var")
    initial_load = zeros(size(force_ratio));
end

static_settings = zeros(1,4);
static_settings(1) = Static_Opts.initial_time_increment;
static_settings(2) = Static_Opts.total_step_time;
static_settings(3) = Static_Opts.minimum_time_increment;
static_settings(4) = Static_Opts.maximum_time_increment;


new_job = JOB_NAME + "_" + job_id(1);

if size(job_id,2) == 2
    new_job = new_job + "_" + job_id(2);
end

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
    end

    if strfind(static_template{iLine,1},'*NODE PRINT,SUMMARY=NO,FREQUENCY = INC_HERE')
        % static_template{iLine,1} = "*NODE PRINT,SUMMARY=NO,FREQUENCY = " + max_inc;
        disp_print_line = iLine;
    end

    if strfind(static_template{iLine,1},'*ENERGY PRINT, FREQUENCY = INC_HERE')
        % static_template{iLine,1} = "*ENERGY PRINT, FREQUENCY = " + max_inc;
        energy_print_line = iLine;
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
        validating_force = (Model.mass*h_evecs).*Static_Opts.perturbation_scale_factor;
end

%-------------------------------------------------------------------------%
geometry = load_geometry(Model);

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
num_nodes = (all_dofs/num_dimensions);
for iNode = 1:num_nodes
    node_label = instance_name + "." + iNode;
    force_label(iNode) = node_label;
    
    for iDimension = 1:num_dimensions
        dimension_label = "," + iDimension + ",";
        force_label(iNode+num_nodes*(iDimension-1),1) = node_label + dimension_label;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
force_transform = Model.mass*Model.reduced_eigenvectors;
coordinate_index = ((1:num_nodes)-1)*num_dimensions;


total_static_steps = sum(num_loadcases);

log_message = sprintf("job " + job_id(1) + ": %i loadcases over %i SEPs" ,[total_static_steps,num_seps]);
logger(log_message,3)

modal_force = zeros(length(Model.reduced_modes),total_static_steps);
sep_id = zeros (1,total_static_steps);
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
            for iDimension = 1:num_dimensions
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
sep_ends = zeros(num_seps,1);


load_step_counter = 0;
total_step_counter = 0;

try
    input_ID = fopen("temp\" + new_job + ".inp","w");
    fprintf(input_ID,'%s\r\n',geometry{:,1});

    for iSep = 1:num_seps
        num_sep_loadcases = num_loadcases(iSep);
        applied_force = force_ratio/num_sep_loadcases;
        force_scale_factors = (1:num_sep_loadcases)';

        if any(restart_type > 0)
            %first "restart" step seems to be problematic
            force_scale_factors = force_scale_factors - 1;
        end


        sep_force = applied_force(:,iSep);
        sep_base_force = initial_load(:,iSep);

        if all(sep_base_force == 0)
            sep_force = sep_force - sep_base_force;
        end
        
        physical_force_bc = force_transform*sep_force;
        physical_force = zeros(all_dofs,1);
        physical_force(Model.node_mapping(:,1),1) = physical_force_bc(Model.node_mapping(:,2),1);

        physical_base_force_bc = force_transform*sep_base_force;
        physical_base_force = zeros(all_dofs,1);
        physical_base_force(Model.node_mapping(:,1),1) = physical_base_force_bc(Model.node_mapping(:,2),1);

        for iLoad = 1:num_sep_loadcases
            load_step_counter = load_step_counter + 1;
            total_step_counter = total_step_counter + 1;
            
            modal_force(:,load_step_counter) = sep_force*force_scale_factors(iLoad) + sep_base_force;
            sep_id(:,load_step_counter) = iSep;

            step_force = physical_force*force_scale_factors(iLoad) + physical_base_force;
           
            
            step_force_label = strings(all_dofs,1);
            for iDimension = 1:num_dimensions
                dimension_span = (1:num_nodes)+(iDimension-1)*num_nodes;
                step_force_label(dimension_span,1) = force_label(dimension_span,1) + step_force(coordinate_index+iDimension,1);
             end

            if iLoad == 1
                step_max_inc = max_inc + num_sep_loadcases;
                load_op = "NEW";
            else
                step_max_inc = max_inc;
                load_op = "MOD";
            end

            static_step = static_template;
            static_step{step_def_line,1} = "*Step, name=STATIC_STEP_" + load_step_counter + ", nlgeom=YES, extrapolation=NO, inc=" + step_max_inc;

            static_step{load_def_line-1} = "*Cload, OP = " + load_op;
            static_step{disp_print_line} = "*NODE PRINT,SUMMARY=NO,FREQUENCY = " + step_max_inc;
            static_step{energy_print_line} = "*ENERGY PRINT, FREQUENCY = " + step_max_inc;

            if deactivated_dofs
                step_force_label(EXCLUDED_DOF) = [];
            end

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
        sep_ends(iSep) = load_step_counter;
    end
catch caught_error
    fclose(input_ID); %ensures file is always closed
    rethrow(caught_error)
end
fclose(input_ID);

setup_time = toc(setup_time_start);
log_message = sprintf("job " + job_id(1) + ": Static input file created: %.1f seconds" ,setup_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_cpus = Model.Static_Options.num_fe_cpus;

abaqus_time_start = tic;

[status,cmdout] = system("cd temp & abaqus job=" + new_job + " cpus=" + num_cpus); %#ok<ASGLU>

while ~isfile("temp\" + new_job + ".dat")
    pause(0.1)
end

while isfile("temp\" + new_job + ".lck")
    pause(0.1)
end

abaqus_time = toc(abaqus_time_start);
log_message = sprintf("job " + job_id(1) + ": Abaqus static analysis complete: %.1f seconds" ,abaqus_time);
logger(log_message,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_processing_time_start = tic;

[displacement,E,additional_data,additional_data_time] = read_abaqus_static_data(new_job,step_type,num_nodes,num_dimensions);

displacement_bc = displacement(Model.node_mapping(:,1),:);
disp_transform = force_transform';

r = disp_transform*displacement_bc;
% theta = displacement_bc - Model.reduced_eigenvectors*r;
theta = displacement_bc;

f = modal_force;

% point = 130;
%  figure;plot([0,r],[0,theta(point,:)],'rx');hold on;plot([0,r],[0,displacement_bc(point,:)],'kx');
step_list = 1:length(E);
final_sep_energy = E(sep_ends);



if clean_data
    remove_index = E > Model.fitting_energy_limit;
    sep_starts = sep_ends - num_loadcases' + 1;
    for iSep = 1:num_seps
        %leave one point ove energy limit in each sep
        remove_sep_index = remove_index(sep_starts(iSep):sep_ends(iSep));
        num_excessive_points = sum(remove_sep_index);
        if num_excessive_points > 0
            remove_sep_index(end+1-num_excessive_points) = 0;
            remove_index(sep_starts(iSep):sep_ends(iSep)) = remove_sep_index;
        end
    end
    
    if all(restart_type ~= 0)
        %first "restart" step seems to be problematic
        remove_index(sep_starts) = 1;
    end
    
    r(:,remove_index) = [];
    theta(:,remove_index) = [];
    E(:,remove_index) = [];
    f(:,remove_index) = [];
    sep_id(:,remove_index) = [];
end

stiff_restarted_seps_index  = (final_sep_energy(restart_type ~= 2) < Model.fitting_energy_limit);

%---
num_points = size(r,2);
point_distance = zeros(1,num_points);
soft_restart_start = [];
soft_restart_end = [];
soft_restarted_seps_index = false(1,num_seps);
for iSep = 1:num_seps
    continue
    if restart_type(iSep) ~= 0 % dont allow multiple softening restarts
        continue
    end
    sep_index = find(sep_id == iSep);
    sep_force = abs(f(1,sep_index));
    [~,sep_order] = sort(sep_force,"ascend");
    sep_index = sep_index(sep_order);
    num_sep_points = length(sep_index);
    previous_point = zeros(size(r,1),1); %NOT CORRECT FOR DOUBLE RESTART
    previous_distance = inf;
    restart_window = 0;
    for iPoint = 1:num_sep_points
        point = r(:,sep_index(iPoint));
        distance = sum((point-previous_point).^2);
        point_distance(sep_index(iPoint)) = distance;
        
        if previous_distance < distance
            if ~restart_window
                restart_window = 1;
                soft_restarted_seps_index(iSep) = true;
                soft_restart_start = [soft_restart_start,sep_index(iPoint-1)];
                end_index = min(iPoint+1,num_sep_points);
                soft_restart_end = [soft_restart_end,sep_index(end_index)];
            else
                end_index = min(iPoint+1,num_sep_points);
                soft_restart_end(end) = sep_index(end_index);
            end
        else
            if restart_window
                restart_window = 0;
            end
        end

        previous_distance = distance;
        previous_point = point;
    end
end
% add points where the point_distance increases dramatically  ->
% significant local softening


%---

data_processing_time = toc(data_processing_time_start) - additional_data_time;
log_message = sprintf("job " + job_id(1) + ": Static data processed: %.1f seconds" ,data_processing_time);
logger(log_message,3)

additional_data_processing_time_start = tic;

switch add_data_type
    case "none"
        additional_data = [];
        
    case "stiffness"
        if clean_data
            step_list(:,remove_index) = [];
        end
        additional_data = parse_stiffness(step_list,new_job,Model.num_dof);

    case "perturbation"
        additional_data = additional_data(Model.node_mapping(:,1),:,:);
        if clean_data
            additional_data(:,:,remove_index) = [];
        end
end
additional_data_processing_time = toc(additional_data_processing_time_start) + additional_data_time;
if add_data_type ~= "none"
    log_message = sprintf("job " + job_id(1) + ": " + add_data_type + " data processed: %.1f seconds" ,additional_data_processing_time);
    logger(log_message,3)
end

%--------------------------
restart_sep = 0;
restarted_loadcases = [];
next_initial_load = [];
next_restart_type = [];
num_restart_loadcases = [];

if any(stiff_restarted_seps_index == 1)
    restart_sep = 1;
    stiff_restarted_loadcases = force_ratio(:,stiff_restarted_seps_index) + initial_load(:,stiff_restarted_seps_index);
    stiff_initial_load = stiff_restarted_loadcases;

    restarted_loadcases = [restarted_loadcases,stiff_restarted_loadcases];
    next_initial_load = [next_initial_load,stiff_initial_load];

    num_stiff_restarts = size(stiff_restarted_loadcases,2);
    next_restart_type = [next_restart_type,ones(1,num_stiff_restarts)];

    num_restart_loadcases = [num_restart_loadcases,num_loadcases(stiff_restarted_seps_index)];
end


if ~isempty(soft_restart_start)
    restart_sep = 1;
    soft_restarted_loadcases = f(:,soft_restart_end);
    soft_initial_load = f(:,soft_restart_start);

    restarted_loadcases = [restarted_loadcases,soft_restarted_loadcases];
    next_initial_load = [next_initial_load,soft_initial_load];

    num_soft_restarts = size(soft_restarted_loadcases,2);
    next_restart_type = [next_restart_type,ones(1,num_soft_restarts)*2];

    num_restart_loadcases = [num_restart_loadcases,num_loadcases(soft_restarted_seps_index)];
end


if restart_sep
    restarted_seps = [find(stiff_restarted_seps_index),find(soft_restarted_seps_index)];
   
    log_message = sprintf("job " + job_id(1) + ": %u/%u SEPs restarted" ,[length(restarted_seps),num_seps]);
    logger(log_message,3)
    
    if size(job_id,2) == 1
        job_id(2) = 2;
    else
        job_id(2) = job_id(2) + 1;
    end
    [r_restart,theta_restart,f_restart,E_restart,additional_data_restart,sep_id_restart] = add_sep_abaqus(restarted_loadcases, ...
        num_restart_loadcases,Static_Opts,max_inc,add_data_type,clean_data,Model,job_id,next_initial_load,next_restart_type);
    
    num_r = size(r,2);
    r_all = [r,r_restart];
    [~,unique_index] = uniquetol(r_all',"ByRows",true);
    unique_index = sort(unique_index,"ascend");
    unique_restart_index = unique_index((num_r+1):end)-num_r;
    r = [r,r_restart(:,unique_restart_index)];

    theta = [theta,theta_restart(:,unique_restart_index)];
    f = [f,f_restart(:,unique_restart_index)];
    E = [E,E_restart(:,unique_restart_index)];
    if ~isempty(additional_data_restart)
        additional_data_restart = additional_data_restart(:,:,unique_restart_index);
    end
    additional_data = cat(3,additional_data,additional_data_restart);
    sep_id = [sep_id,restarted_seps(sep_id_restart(:,unique_restart_index))];


end

end

