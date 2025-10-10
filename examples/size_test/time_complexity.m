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
seed_sizes = [0.00307,0.0023,0.002,0.00174,0.00163,0.00157,0.001481,0.00142,0.00138,0.00136,0.00133,0.001293,0.00129];
num_workers = 4;
num_static_repeats = 2;

num_dynamic_repeats = 10;
%-------
data_path = "data\size_data";

static_log_lines = [
    "mems_arch:", "Eigenvectors:";
    "Eigenvectors:","Model Initialised:";
    "Model Initialised:","Dataset scaffold created:";
    "Dataset scaffold created:", "Verification step 1";
    "Verification step 1", "Verification step 2"
    ];

dynamic_log_lines = @(id) [
    "Dynamic analysis " + id, "Dynamic Dataset initalised:";
    "Dynamic Dataset initalised:", "EoM precomputations:";
    "EoM precomputations:", "Linear solution found";
    "Linear solution found:", "EoM loaded:";
    "EoM loaded:", "Backbone:";
    "Backbone:","EoM loaded:";
    "EoM loaded:", "Initial condition sweep:";
    "Initial condition sweep:", "EoM loaded:";
    "EoM loaded:", "Backbone:"
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
Dynamic_Data.seed_sizes = seed_sizes;

num_seeds = length(seed_sizes);
num_dof = zeros(1,num_seeds);

total_time = zeros(1,num_seeds);
matrix_time = zeros(1,num_seeds);
initialisation_time = zeros(1,num_seeds);
scaffold_time = zeros(1,num_seeds);
verification_time = zeros(2,num_seeds);
free_static_memory = cell(1,num_seeds);

free_dynamic_memory = cell(1,num_seeds);
dynamic_time = zeros(3,num_seeds);

for iSeed = 1:num_seeds
    seed_size = seed_sizes(iSeed);

    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);

    for iRepeat = 1:num_static_repeats
        create_parallel_pool(0);
        create_parallel_pool(num_workers);
        clear Static_Data
        clear Model
        clear Dyn_Data

        delete_static_data(get_system_name(system_name,modes));
        delete_cache(system_name,"force",energy_limit)
        delete_cache(system_name,"matrices")
        start_memory_profiler("display",0)
        %create static data
        total_time_start = tic;
        Model = Dynamic_System(system_name,energy_limit,modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);
        Static_Data = Static_Dataset(Model);
        Static_Data.save_data;
        total_time(1,iSeed) = toc(total_time_start);

        stop_memory_profiler
        [memory_data,memory_duration] = get_free_memory;
        Memory.data = memory_data;
        Memory.duration = memory_duration;
        free_static_memory{1,iSeed} = Memory;


        log_data = read_log(static_log_lines);
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
        Size_Data(iRepeat).free_memory(1,iSeed) = free_static_memory(1,iSeed);
        Size_Data(iRepeat).num_dofs(1,iSeed) = num_dof(1,iSeed);

        save(data_path + "\size_data","Size_Data")
        Model.save_log
    end

    for iRepeat = 1:num_dynamic_repeats
        clear Static_Data
        clear Model
        clear Dyn_Data
        delete_dynamic_data(get_system_name(system_name,modes));

        start_memory_profiler("display",0)

   
        log_message = sprintf("Dynamic analysis %u/%u" ,iRepeat,num_dynamic_repeats);
        logger(log_message,1)

        dynamic_time_one_start = tic;
        Dyn_Data = initalise_dynamic_data(get_system_name(system_name,modes));
        %--
        % Additional_Output.output = "physical displacement";
        % Additional_Output.type = "max";
        % Additional_Output.dof = 66539;
        % Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
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
        potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2.69e6,[1.5e-7,1e-7],"figure",0);
        dynamic_time(2,iSeed) = toc(dynamic_time_two_start);
        %---
        dynamic_time_three_start = tic;
        Continuation_Opts.collation_degree = 10;
        Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
        dynamic_time(3,iSeed) = toc(dynamic_time_three_start);
        %---
        stop_memory_profiler
        [memory_data,memory_duration] = get_free_memory;
        Memory.data = memory_data;
        Memory.duration = memory_duration;
        free_dynamic_memory{1,iSeed} = Memory;
        %--
        Dynamic_Data(iRepeat).dynamic_time(:,iSeed) = dynamic_time(:,iSeed);
        Dynamic_Data(iRepeat).free_memory(1,iSeed) = free_dynamic_memory(1,iSeed);
        Dynamic_Data(iRepeat).num_dofs(1,iSeed) = num_dof(1,iSeed);
        
        %-----------------------------
        log_data = read_log(dynamic_log_lines(iRepeat));

        eom_construction_time = sum(log_data([1,2,4,6,8]));
        backbone_1_time = sum(log_data([3,5]));
        ic_time = log_data(7);
        backbone_2_time = log_data(9);

        Dynamic_Data(iRepeat).eom_construction_time(:,iSeed) = eom_construction_time;
        Dynamic_Data(iRepeat).backbone_1_time(:,iSeed) = backbone_1_time;
        Dynamic_Data(iRepeat).ic_time(:,iSeed) = ic_time;
        Dynamic_Data(iRepeat).backbone_2_time(:,iSeed) = backbone_2_time;
        %-----------------------------

        save(data_path + "\dynamic_data","Dynamic_Data")

        
        % Model.save_log
    end
end
