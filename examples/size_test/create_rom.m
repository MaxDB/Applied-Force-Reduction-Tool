clear
close all

num_iterations = 1;

%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 0.8;
initial_modes = [1];
added_modes = [6,11];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.output_format = "binary";
Static_Opts.num_loadcases = 5;
Static_Opts.num_validation_modes = 20;
%------------------------------------------%
Calibration_Opts.Static_Opts = struct([]);
Calibration_Opts.calibration_scale_factor = 1.5;
Verification_Opts.num_added_points = 5;
Verification_Opts.maximum_iterations = 1;

if isempty(gcp('nocreate')) && Static_Opts.max_parallel_jobs > 1
    parpool;
end

rom_one_base = zeros(1,num_iterations);
rom_one_validation_data = zeros(1,num_iterations);
rom_one_orbits = zeros(1,num_iterations);
rom_one_orbit_validation = zeros(1,num_iterations);

rom_two_base = zeros(1,num_iterations);
rom_two_validation = zeros(1,num_iterations);
rom_two_orbits = zeros(1,num_iterations);
rom_two_orbit_validation = zeros(1,num_iterations);

for iCount = 1:num_iterations
    close all
    delete_static_data(system_name+"_"+initial_modes);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    % base_time_start = tic;
    % Static_Opts.additional_data = "none";
    % Static_Data_Base = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts);
    % rom_one_base(iCount) = toc(base_time_start);

    % base_time_start = tic;
    % Static_Data_Base_2 = two_mode_rom(Static_Data_Base,added_modes);
    % rom_two_base(iCount) = toc(base_time_start);

    clear("Static_Data_Base")
    clear("Static_Data_Base_2")
    delete_static_data(system_name+"_"+initial_modes);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    validation_time_start = tic;
    Static_Opts.additional_data = "stiffness";
    Static_Data_Validation = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts);
    rom_one_validation_data(iCount) = toc(validation_time_start);

    % orbit_time_start = tic;
    % Dyn_Data = one_mode_rom_orbits(system_name+"_"+initial_modes);
    % rom_one_orbits(1,iCount) = toc(orbit_time_start);
    % 
    % orbit_validation_start = tic;
    % Dyn_Data = one_mode_rom_validation(Dyn_Data);
    % rom_one_orbit_validation(1,iCount) = toc(orbit_validation_start);


    % clear("Dyn_Data")
    % close all
    % %----------------------------------------------------------------
    validation_time_start = tic;
    Static_Data_Validation = two_mode_rom(Static_Data_Validation,added_modes(1));
    rom_two_validation(1,iCount) = toc(validation_time_start);
    % 
    % orbit_time_start = tic;
    % Dyn_Data = two_mode_rom_orbits(system_name+"_"+initial_modes + added_modes);
    % rom_two_orbits(1,iCount) = toc(orbit_time_start);
    % 
    % orbit_validation_start = tic;
    % Dyn_Data = two_mode_rom_validation(Dyn_Data);
    % rom_two_orbit_validation(1,iCount) = toc(orbit_validation_start);
    % 
    % %----------------------------------------------------------------
    % clear("Dyn_Data")
    % clear("Static_Data_Validation")
   
end

print_mean_time(rom_one_base,"Static Data")
print_mean_time(rom_one_validation_data,"Validation Data")
print_mean_time(rom_one_validation_data-rom_one_base,"Data diff")
print_mean_time(rom_one_orbits,"Orbits")
print_mean_time(rom_one_orbit_validation,"Orbit validation")

total_one = rom_one_validation_data + sum(rom_one_orbits,1) + sum(rom_one_orbit_validation,1);
print_mean_time(total_one,"Total")
fprintf("---\n\n");

print_mean_time(rom_two_base,"Static Data")
print_mean_time(rom_two_validation,"Validation Data")
print_mean_time(rom_two_validation-rom_two_base,"Data diff")
print_mean_time(rom_two_orbits,"Orbits")
print_mean_time(rom_two_orbit_validation,"Orbit validation")

total_two = rom_two_validation + rom_two_orbits + rom_two_orbit_validation;
print_mean_time(total_two,"Total")
fprintf("---\n\n");

%-----------------------------------
function Static_Data = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts)


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts,"calibration_opts",Calibration_Opts);
Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;
end

function Dyn_Data = one_mode_rom_orbits(system_name)
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = "all";
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% 
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e0;
Continuation_Opts.max_inc = 1e0;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
end

function Dyn_Data = one_mode_rom_validation(Dyn_Data)
compare_validation(Dyn_Data,"validation error",1,"all");
end


function Static_Data = two_mode_rom(Static_Data,added_modes)
Static_Data = Static_Data.update_model(added_modes);
Static_Data = Static_Data.create_dataset;
Static_Data.save_data;
end

function Dyn_Data = two_mode_rom_orbits(system_name)
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 5e-2;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 8;
%-----------------------------------------%
Dyn_Data = initalise_dynamic_data(system_name);

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
end

function Dyn_Data = two_mode_rom_validation(Dyn_Data)
compare_validation(Dyn_Data,"validation error",1,1:10);
end