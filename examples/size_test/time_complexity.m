clear
close all
%---------------------------------
% 0.01    --> 7,000
% 0.005   --> 24,500
% 0.0044 is too big
% 0.003   --> 100,000
% 0.002   --> 300,000
% 0.00125 --> 1,000,000
% 0.001   --> 2,000,000
% seed_sizes = [0.0044,0.003,0.00195,0.00125];

% dof ≈ 0.0458 * seed_size ^ -2.53
%-------------------------------
% seed_sizes = [0.003,0.0024,0.0021,0.0018,0.0017,0.0016,0.0015,0.0014,0.0013,0.00125];
seed_sizes = [0.003,0.0024];
num_workers = 4;
num_repeats = 2;

dynamic_data = 1;
%-------
data_path = "data\size_data";

log_lines = [
    "mems_arch:", "Eigenvectors:";
    "Eigenvectors:","Model Initialised:";
    "Model Initialised:","Dataset scaffold created:";
    "Dataset scaffold created:", "Verification step 1";
    "Verification step 1", "Verification step 2"
    ];
%-----
set_logging_level(1)
set_visualisation_level(0)
system_name = "mems_arch";
energy_limit = 0.8;
modes = [1,6];

Calibration_Opts.calibration_scale_factor = 1;
%----
Static_Opts.max_parallel_jobs = num_workers;
create_parallel_pool(num_workers);

data_path = data_path + "\workers_" + num_workers;

%-- setup
if ~isfolder(data_path)
    mkdir(data_path)
end
%--
Size_Data.seed_sizes = seed_sizes;

num_seeds = length(seed_sizes);
num_dof = zeros(1,num_seeds);

total_time = zeros(1,num_seeds);
matrix_time = zeros(1,num_seeds);
initialisation_time = zeros(1,num_seeds);
scaffold_time = zeros(1,num_seeds);
verification_time = zeros(2,num_seeds);
free_memory = cell(1,num_seeds);

dynamic_time = zeros(3,num_seeds);

for iSeed = 1:num_seeds
    seed_size = seed_sizes(iSeed);

    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);

    for iRepeat = 1:num_repeats
        clear Static_Data
        clear Model
        clear Dyn_Data

        delete_static_data(get_system_name(system_name,modes));
        delete_cache(system_name,"force",energy_limit)
        delete_cache(system_name,"matrices")
        start_memory_profiler
        %create static data
        total_time_start = tic;
        Model = Dynamic_System(system_name,energy_limit,modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);
        Static_Data = Static_Dataset(Model);
        Static_Data.save_data;
        total_time(1,iSeed) = toc(total_time_start);

        stop_memory_profiler
        memory_data = get_free_memory;
        free_memory{1,iSeed} = memory_data;

        log_data = read_log(log_lines);
        matrix_time(1,iSeed) = log_data(1);
        initialisation_time(1,iSeed) = log_data(2);
        scaffold_time(1,iSeed) = log_data(3);
        verification_time(1,iSeed) = log_data(4);
        verification_time(2,iSeed) = log_data(5);

        Size_Data(iRepeat).total_time(1,iSeed) = total_time(1,iSeed); %#ok<*SAGROW>
        Size_Data(iRepeat).matrix_time(1,iSeed) = matrix_time(1,iSeed);
        Size_Data(iRepeat).initialisation_time(1,iSeed) = initialisation_time(1,iSeed);
        Size_Data(iRepeat).scaffold_time(1,iSeed) = scaffold_time(1,iSeed);
        Size_Data(iRepeat).verification_time(:,iSeed) = verification_time(:,iSeed);
        Size_Data(iRepeat).free_memory(1,iSeed) = free_memory(1,iSeed);
        Size_Data(iRepeat).num_dofs(1,iSeed) = num_dof(1,iSeed);

        save(data_path + "\size_data","Size_Data")
        
        if ~dynamic_data
            continue
        end
        
        dynamic_time_one_start = tic;
        Dyn_Data = initalise_dynamic_data(get_system_name(system_name,modes));
        %--
        Additional_Output.output = "physical displacement";
        Additional_Output.type = "max";
        Additional_Output.dof = 66539;
        Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
        % --------- Continuation Settings ---------%
        Continuation_Opts.initial_inc = 5e-1;
        Continuation_Opts.max_inc = 5e-1;
        Continuation_Opts.min_inc = 1e-2;
        Continuation_Opts.forward_steps = 100;
        Continuation_Opts.backward_steps = 0;
        Continuation_Opts.collation_degree = 8;
        % ----------------------------------------- %
        Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
        dynamic_time(1,iSeed) = toc(dynamic_time_one_start);
        %---
        dynamic_time_two_start = tic;
        potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2.69e6,[1.5e-7,1e-7]);
        dynamic_time(2,iSeed) = toc(dynamic_time_two_start);
        %---
        dynamic_time_three_start = tic;
        Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
        dynamic_time(3,iSeed) = toc(dynamic_time_three_start);
        %---

        Size_Data(iRepeat).dynamic_time(:,iSeed) = dynamic_time(:,iSeed);
        save(data_path + "\size_data","Size_Data")

    end
end
