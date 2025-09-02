clear
close all

num_iterations = 1;

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

initial_time = zeros(1,num_iterations);
calibration = zeros(3,num_iterations);
initial_points = zeros(3,num_iterations);
verification_step_1 = zeros(3,num_iterations);
verification_step_2 = zeros(3,num_iterations);
verification_step_3 = zeros(3,num_iterations);

for iCount = 1:num_iterations
    close all
    delete_static_data(system_name+"_"+initial_modes);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    initial_time_start = tic;
    create_model(system_name,0,initial_modes,Static_Opts);
    initial_time(iCount) = toc(initial_time_start);

    calibration_time_start = tic;
    Model = create_model(system_name,energy_limit,initial_modes,Static_Opts);
    calibration(1,iCount) = toc(calibration_time_start);

    Static_Data = Static_Dataset(Model,Verification_Opts);

   clear Model
   clear Static_Dataset
end

print_mean_time(rom_one_base,"Static Data")
print_mean_time(rom_one_validation_data,"Validation Data")
print_mean_time(rom_one_validation_data-rom_one_base,"Data diff")
print_mean_time(rom_one_orbits(1,:),"Orbits 1")
print_mean_time(rom_one_orbit_validation(1,:),"Orbits 1 validation")
print_mean_time(rom_one_orbits(2,:),"Orbits 2")
print_mean_time(rom_one_orbit_validation(2,:),"Orbits 2 validation")

total_one = rom_one_validation_data + sum(rom_one_orbits,1) + sum(rom_one_orbit_validation,1);
print_mean_time(total_one,"Total")
fprintf("---\n\n");



%-----------------------------------

