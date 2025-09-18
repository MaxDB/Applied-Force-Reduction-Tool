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
% seed_sizes = [0.003,0.00195,0.00125];#
seed_sizes = 0.003;
num_workers = 4;


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

for iSeed = 1:num_seeds
    seed_size = seed_sizes(iSeed);

    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);
    
    start_memory_profiler
    %create static data
    total_time_start = tic;
    Model = Dynamic_System(system_name,energy_limit,modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);
    Static_Data = Static_Dataset(Model);
    Static_Data.save_data;
    total_time(1,iSeed) = toc(total_time_start = tic);

    stop_memory_profiler
    memory_data = get_free_memory;
    free_memory{1,iSeed} = memory_data;

    log_data = read_log(log_lines);
    matrix_time(1,iSeed) = log_data(1);
    initialisation_time(1,iSeed) = log_data(2);
    scaffold_time(1,iSeed) = log_data(3);
    verification_time(1,iSeed) = log_data(4);
    verification_time(2,iSeed) = log_data(5);

    Size_Data.total_time = total_time;
    Size_Data.matrix_time = matrix_time;
    Size_Data.initialisation_time = initialisation_time;
    Size_Data.scaffold_time = scaffold_time;
    Size_Data.verification_time = verification_time;
    Size_Data.free_memory = free_memory;

    save(data_path + "\size_data","Size_Data") 
end
