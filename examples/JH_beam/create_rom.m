clear
close all

num_iterations = 3;

%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "JH_beam_2d";
energy_limit = 0.015; 
initial_modes = 1;
added_modes = [3,5];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
%------------------------------------------%
if isempty(gcp('nocreate')) && Static_Opts.max_parallel_jobs > 1
    parpool;
end

log_lines = [
    "Eigenvectors:", "Model Initialised:";
    "Model Initialised:","Dataset scaffold created:";
    "Dataset scaffold created:","Verification step 1";
    "Verification step 1", "Verification step 2";
    "Verification step 2", "Verification step 3"
    ];


rom_one = get_system_name(system_name,initial_modes);
rom_two = get_system_name(system_name,[initial_modes,added_modes(1)]);
rom_three = get_system_name(system_name,[initial_modes,added_modes]);

initial_time = zeros(1,num_iterations);
calibration = zeros(3,num_iterations);
initial_points = zeros(3,num_iterations);
verification_step_1 = zeros(3,num_iterations);
verification_step_2 = zeros(3,num_iterations);
verification_step_3 = zeros(3,num_iterations);

total_time = zeros(3,num_iterations);

for iCount = 1:num_iterations
    close all
   
    delete_static_data(rom_one);
    delete_static_data(rom_two);
    delete_static_data(rom_three);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")
    
    total_time_start = tic;
    initial_time_start = tic;
    Dynamic_System(system_name,0,initial_modes,"static_opts",Static_Opts);
    initial_time(iCount) = toc(initial_time_start);

    calibration_time_start = tic;
    Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
    calibration(1,iCount) = toc(calibration_time_start);

    Static_Data = Static_Dataset(Model);
    Static_Data.save_data;
    total_time(1,iCount) = toc(total_time_start);

    log_data = read_log(log_lines);
    calibration(1,iCount) = log_data(1);
    initial_points(1,iCount) = log_data(2);
    verification_step_1(1,iCount) = log_data(3);


    %-------------------------------------------
    total_time_start = tic;
    Static_Data = Static_Data.update_model(added_modes(1));
    Static_Data = Static_Data.create_dataset;
    Static_Data.save_data;
    total_time(2,iCount) = toc(total_time_start);

    log_data = read_log(log_lines);
    calibration(2,iCount) = log_data(1);
    initial_points(2,iCount) = log_data(2);
    verification_step_1(2,iCount) = log_data(3);
    verification_step_2(2,iCount) = log_data(4);
    verification_step_3(2,iCount) = log_data(5);
    %-------------------------------------------
    total_time_start = tic;
    
    
    calibration_time_start = tic;
    Model = Dynamic_System(system_name,energy_limit,[initial_modes,added_modes],"static_opts",Static_Opts);
    calibration(3,iCount) = toc(calibration_time_start);
    Static_Data = Static_Dataset(Model);
    Static_Data.save_data;
    total_time(3,iCount) = toc(total_time_start);

    log_data = read_log(log_lines);
    calibration(3,iCount) = log_data(1);
    initial_points(3,iCount) = log_data(2);
    verification_step_1(3,iCount) = log_data(3);
    verification_step_2(3,iCount) = log_data(4);
    %-------------------------------------------
    clear Model
    clear Static_Dataset
end

print_mean_time(initial_time,"Initialisation")
fprintf("\n");
print_mean_time(calibration,"Calibration")
fprintf("\n");
print_mean_time(initial_points ,"Initial points")
fprintf("\n");
print_mean_time(verification_step_1,"Verification iteration 1")
fprintf("\n");
print_mean_time(verification_step_2,"Verification iteration 2")
fprintf("\n");
print_mean_time(verification_step_3,"Verification iteration 3")
fprintf("\n");
print_mean_time(total_time,"Total ROM creation")
% 
% total_one = rom_one_validation_data + sum(rom_one_orbits,1) + sum(rom_one_orbit_validation,1);
% print_mean_time(total_one,"Total")



%-----------------------------------

