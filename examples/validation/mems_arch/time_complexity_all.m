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
% seed_sizes = [0.00307,0.0023,0.002,0.00174,0.00163,0.00157,0.001481,0.00142,0.00138,0.00136,0.00133,0.001293,0.00129];
seed_sizes = [0.00133];

num_workers = 2;
num_repeats = 2;

%-------
data_path = "data\size_data";

static_log_lines_one = [
    "mems_arch:", "Eigenvectors:";
    "Eigenvectors:","Model Initialised:";
    "Model Initialised:","Dataset scaffold created:";
    "Dataset scaffold created:", "Verification step 1";
    "Verification step 1", "ROM Dataset Created"
    ];

static_log_lines_two = [
    "mems_arch:", "Eigenvectors:";
    "Eigenvectors:","Model Initialised:";
    "Model Initialised:","Dataset scaffold created:";
    "Dataset scaffold created:", "Verification step 1";
    "Verification step 1", "Verification step 2";
    "Verification step 2", "ROM Dataset Created"
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
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 20;
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

total_time = zeros(2,num_seeds);
matrix_time = zeros(2,num_seeds);
initialisation_time = zeros(2,num_seeds);
scaffold_time = zeros(2,num_seeds);
verification_time = zeros(2,num_seeds);
perturbation_time = zeros(2,num_seeds);
free_static_memory = cell(1,num_seeds);

validation_time = zeros(2,num_seeds);

free_dynamic_memory = cell(1,num_seeds);
dynamic_time = zeros(3,num_seeds);

% start_crash_monitoring()

for iSeed = 1:num_seeds
    seed_time_start = tic;
    seed_size = seed_sizes(iSeed);

    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);

    for iRepeat = 1:num_repeats
        create_parallel_pool(0);

        clear Static_Data
        clear Model
        clear Dyn_Data
        close all

        try
            create_parallel_pool(num_workers);
        catch
            pause(30)
            create_parallel_pool(num_workers);
        end
        

        delete_static_data(get_system_name(system_name,modes(1)));
        delete_static_data(get_system_name(system_name,modes));
        delete_cache(system_name,"force",energy_limit)
        delete_cache(system_name,"matrices")
        start_memory_profiler("display",0)
        %&&create static data
        total_time_start = tic;
        % one mode
        Model = Dynamic_System(system_name,energy_limit,modes(1),"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);
        Static_Data = Static_Dataset(Model);
        Static_Data.save_data;
        total_time(1,iSeed) = toc(total_time_start);

        log_data = read_log(static_log_lines_one);
        matrix_time(1,iSeed) = log_data(1);
        initialisation_time(1,iSeed) = log_data(2);
        scaffold_time(1,iSeed) = log_data(3);
        verification_time(1,iSeed) = log_data(4);
        perturbation_time(1,iSeed) = log_data(5);
        %------
        validation_time_start = tic;
        one_mode_validation()
        validation_time(1,iSeed) = toc(validation_time_start);
        %------------------
        total_time_start = tic;
        %two mode
        Static_Data = Static_Data.update_model(modes(2));
        Static_Data = Static_Data.create_dataset;
        Static_Data.save_data;
        total_time(2,iSeed) = toc(total_time_start);
        
        %------
        validation_time_start = tic;
        two_mode_validation()
        validation_time(2,iSeed) = toc(validation_time_start);
        %------------------
        

        stop_memory_profiler
        [memory_data,memory_duration] = get_free_memory;
        Memory.data = memory_data;
        Memory.duration = memory_duration;
        free_static_memory{1,iSeed} = Memory;
        
        
        log_data = read_log(static_log_lines_two);
        matrix_time(2,iSeed) = log_data(1);
        initialisation_time(2,iSeed) = log_data(2);
        scaffold_time(2,iSeed) = log_data(3);
        verification_time(2,iSeed) = sum(log_data(4:5));
        perturbation_time(2,iSeed) = log_data(6);

        Size_Data(iRepeat).total_time(:,iSeed) = total_time(:,iSeed); %#ok<*SAGROW>
        Size_Data(iRepeat).matrix_time(:,iSeed) = matrix_time(:,iSeed);
        Size_Data(iRepeat).initialisation_time(:,iSeed) = initialisation_time(:,iSeed);
        Size_Data(iRepeat).scaffold_time(:,iSeed) = scaffold_time(:,iSeed);
        Size_Data(iRepeat).verification_time(:,iSeed) = verification_time(:,iSeed);
        Size_Data(iRepeat).perturbation_time(:,iSeed) = perturbation_time(:,iSeed);
        Size_Data(iRepeat).validation_time(:,iSeed) = validation_time(:,iSeed);
        Size_Data(iRepeat).free_memory(1,iSeed) = free_static_memory(1,iSeed);
        Size_Data(iRepeat).num_dofs(1,iSeed) = num_dof(1,iSeed);

        save(data_path + "\size_data","Size_Data")
        % Model.save_log
    end
    seed_time = toc(seed_time_start);
    simulation_update("MEMS arch validation (%i) in %.1f: %i/%i",[num_dof(iSeed),seed_time/60,iSeed,num_seeds])

end


function one_mode_validation
system_name = "mems_arch_1";
Dyn_Data = initalise_dynamic_data(system_name);
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [0,36e-3,10e-3];
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 5e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 20;
% -----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

compare_validation(Dyn_Data,"validation error",1,1:6);


end

function two_mode_validation
% create_parallel_pool(8);
system_name = "mems_arch_16";
Dyn_Data = initalise_dynamic_data(system_name);
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [0,36e-3,10e-3];
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 5e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 40;
% -----------------------------------------%
try
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
catch
    Continuation_Opts.collocation_degree = 8;
    Continuation_Opts.initial_discretisation_num = 20;
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    Continuation_Opts.collocation_degree = 6;
    Continuation_Opts.initial_discretisation_num = 40;
end

Dyn_Data_One_Mode = initalise_dynamic_data("mems_arch_1");
% Dyn_Data_One_Mode = Dyn_Data_One_Mode.validate_solution(1,6);
Sol = Dyn_Data_One_Mode.load_solution(1,"validation");
unstable_index = find(Sol.h_stability>1.01,3);
restart_index = unstable_index(2);
[orbit,validation_orbit] = Dyn_Data_One_Mode.get_orbit(1,restart_index,1);
[q,q_dot] = Dyn_Data_One_Mode.get_modal_validation_orbit(1,restart_index);
[min_ke,min_index] = min(sum(q_dot.^2,1)); 
test_ic = q(:,min_index);

potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2*pi/orbit.T,test_ic);

Continuation_Opts.backward_steps = 5;
try
    Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
catch
    try
        Continuation_Opts.collocation_degree = 8;
        Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
    catch
        Continuation_Opts.collocation_degree = 15;
        Continuation_Opts.initial_discretisation_num = 100;
        Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
    end
end


compare_validation(Dyn_Data,"validation error",[1,2],"all")

Dyn_Data.validate_solution(1,[5,11,13]);
Dyn_Data.validate_solution(2,[5,11,13],"load_data",1);
end