clear
close all
%---
num_seed_sizes = 100;
max_seed_size = 0.01;
min_seed_size = 0.000035;



%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "clamped_beam";
energy_limit = 0.01;
initial_modes = [1];
added_modes = [3];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.output_format = "binary";
Static_Opts.num_loadcases = 6;
%------------------------------------------%

if isempty(gcp('nocreate')) && Static_Opts.max_parallel_jobs > 1
    parpool;
end

rom_one_base = zeros(1,num_seed_sizes);
rom_one_validation_data = zeros(1,num_seed_sizes);
rom_one_orbits = zeros(2,num_seed_sizes);
rom_one_orbit_validation = zeros(2,num_seed_sizes);

rom_two_base = zeros(1,num_seed_sizes);
rom_two_validation = zeros(1,num_seed_sizes);
rom_two_orbits = zeros(1,num_seed_sizes);
rom_two_orbit_validation = zeros(1,num_seed_sizes);

seed_sizes = logspace(log10(max_seed_size),log10(min_seed_size),num_seed_sizes);
mesh_dof = zeros(1,num_seed_sizes);
for iStep = 1:num_seed_sizes
    close all

    mesh_dof(iStep) = set_seed_size(system_name,seed_sizes(iStep));
   

    base_time_start = tic;
    Static_Opts.additional_data = "none";
    Static_Data_Base = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts);
    rom_one_base(iStep) = toc(base_time_start);

    base_time_start = tic;
    Static_Data_Base_2 = two_mode_rom(Static_Data_Base,added_modes);
    rom_two_base(iStep) = toc(base_time_start);

    clear("Static_Data_Base")
    clear("Static_Data_Base_2")
    delete_static_data(system_name+"_"+initial_modes);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    validation_time_start = tic;
    Static_Opts.additional_data = "stiffness";
    Static_Data_Validation = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts);
    rom_one_validation_data(iStep) = toc(validation_time_start);

    orbit_time_start = tic;
    Dyn_Data = one_mode_rom_orbits(system_name+"_"+initial_modes,1);
    rom_one_orbits(1,iStep) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = one_mode_rom_validation(Dyn_Data,1);
    rom_one_orbit_validation(1,iStep) = toc(orbit_validation_start);

    orbit_time_start = tic;
    Dyn_Data = one_mode_rom_orbits(Dyn_Data,2);
    rom_one_orbits(2,iStep) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = one_mode_rom_validation(Dyn_Data,2);
    rom_one_orbit_validation(2,iStep) = toc(orbit_validation_start);

    clear("Dyn_Data")
    close all
    %----------------------------------------------------------------
    validation_time_start = tic;
    Static_Data_Validation = two_mode_rom(Static_Data_Validation,added_modes);
    rom_two_validation(1,iStep) = toc(validation_time_start);

    orbit_time_start = tic;
    Dyn_Data = two_mode_rom_orbits(system_name+"_"+initial_modes + added_modes);
    rom_two_orbits(1,iStep) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = two_mode_rom_validation(Dyn_Data);
    rom_two_orbit_validation(1,iStep) = toc(orbit_validation_start);
    
    %----------------------------------------------------------------
    clear("Dyn_Data")
    clear("Static_Data_Validation")

     disp(iStep + "/" + num_seed_sizes)
end