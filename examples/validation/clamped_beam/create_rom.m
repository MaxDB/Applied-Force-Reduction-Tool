clear
close all

num_iterations = 1;

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
%------------------------------------------%
create_parallel_pool(Static_Opts.max_parallel_jobs);


rom_one_static_data = zeros(1,num_iterations);
rom_one_orbits = zeros(2,num_iterations);
rom_one_orbit_validation = zeros(2,num_iterations);

rom_two_static_data = zeros(1,num_iterations);
rom_two_orbits = zeros(1,num_iterations);
rom_two_orbit_validation = zeros(1,num_iterations);

for iCount = 1:num_iterations
    close all
    delete_static_data(system_name+"_"+initial_modes);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    delete_static_data(system_name+"_"+initial_modes);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    verification_time_start = tic;
    Static_Opts.additional_data = "stiffness";
    Static_Data = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts);
    rom_one_static_data(iCount) = toc(verification_time_start);

    orbit_time_start = tic;
    Dyn_Data = one_mode_rom_orbits(system_name+"_"+initial_modes,1);
    rom_one_orbits(1,iCount) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = one_mode_rom_validation(Dyn_Data,1);
    rom_one_orbit_validation(1,iCount) = toc(orbit_validation_start);

    orbit_time_start = tic;
    Dyn_Data = one_mode_rom_orbits(Dyn_Data,2);
    rom_one_orbits(2,iCount) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = one_mode_rom_validation(Dyn_Data,2);
    rom_one_orbit_validation(2,iCount) = toc(orbit_validation_start);

    clear("Dyn_Data")
    % close all
    %----------------------------------------------------------------
    verification_time_start = tic;
    Static_Data = two_mode_rom(Static_Data,added_modes);
    rom_two_static_data(1,iCount) = toc(verification_time_start);

    orbit_time_start = tic;
    Dyn_Data = two_mode_rom_orbits(system_name+"_"+initial_modes + added_modes);
    rom_two_orbits(1,iCount) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = two_mode_rom_validation(Dyn_Data);
    rom_two_orbit_validation(1,iCount) = toc(orbit_validation_start);
    
    %----------------------------------------------------------------
    clear("Dyn_Data")
    clear("Static_Data")
   
end

print_mean_time(rom_one_static_data,"Static Data")
print_mean_time(rom_one_orbits(1,:),"Orbits 1")
print_mean_time(rom_one_orbit_validation(1,:),"Orbits 1 validation")
print_mean_time(rom_one_orbits(2,:),"Orbits 2")
print_mean_time(rom_one_orbit_validation(2,:),"Orbits 2 validation")

total_one = rom_one_static_data + sum(rom_one_orbits,1) + sum(rom_one_orbit_validation,1);
print_mean_time(total_one,"Total")
fprintf("---\n\n");

print_mean_time(rom_two_static_data,"Static Data")
print_mean_time(rom_two_orbits,"Orbits")
print_mean_time(rom_two_orbit_validation,"Orbit validation")

total_two = rom_two_static_data + rom_two_orbits + rom_two_orbit_validation;
print_mean_time(total_two,"Total")
fprintf("---\n\n");

%-----------------------------------

