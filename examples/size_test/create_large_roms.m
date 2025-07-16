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
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.output_format = "binary";
Static_Opts.num_loadcases = 3;
Static_Opts.num_validation_modes = 20;
%------------------------------------------%
Calibration_Opts.Static_Opts = struct([]);
Calibration_Opts.calibration_scale_factor = 1.5;
Verification_Opts.num_added_points = 1;
Verification_Opts.maximum_iterations = 1;
Verification_Opts.maximum_interpolation_error =  [1e-2;1e-2];

if isempty(gcp('nocreate')) && Static_Opts.max_parallel_jobs > 1
    parpool;
end


rom_three.name = "mems_arch_1611";
rom_three.modes =  [1,6,11];

rom_three.base = zeros(1,num_iterations);
rom_three.validation_data = zeros(1,num_iterations);
rom_three.orbits = zeros(1,num_iterations);
rom_three.orbit_validation = zeros(1,num_iterations);

rom_four.name = "mems_arch_15611";
rom_four.modes = [1,5,6,11];

rom_four.base = zeros(1,num_iterations);
rom_four.validation = zeros(1,num_iterations);
rom_four.orbits = zeros(4,num_iterations);
rom_four.orbit_validation = zeros(3,num_iterations);


%times assume mass matrix and calibration for 1st and 6th mode has already
%been done
for iCount = 1:num_iterations
    close all
    delete_static_data(rom_three.name);
    delete_static_data(rom_four.name);
    delete_cache(system_name,"force",energy_limit,11) 
    delete_cache(system_name,"force",energy_limit,5) 
    
    base_time_start = tic;
    Static_Opts.additional_data = "none";
    Static_Data_Base = three_mode_rom(system_name,energy_limit,rom_three.modes,Static_Opts,Verification_Opts,Calibration_Opts);
    rom_three.base(iCount) = toc(base_time_start);   

    base_time_start = tic;
    Static_Opts.additional_data = "none";
    Static_Data_Base = four_mode_rom(system_name,energy_limit,rom_four.modes,Static_Opts,Verification_Opts,Calibration_Opts);
    rom_four.base(iCount) = toc(base_time_start);   

    clear("Static_Data_Base")
    delete_static_data(rom_three.name);
    delete_static_data(rom_four.name);
    delete_cache(system_name,"force",energy_limit,11)
    delete_cache(system_name,"force",energy_limit,5) 

    validation_time_start = tic;
    Static_Opts.additional_data = "stiffness";
    Static_Data_Validation = three_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts);
    rom_three.validation_data(iCount) = toc(validation_time_start);

    validation_time_start = tic;
    Static_Opts.additional_data = "stiffness";
    Static_Data_Validation = four_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts);
    rom_four.validation_data(iCount) = toc(validation_time_start);

end
exclude_data = [];
fprintf("---\n\n");
fprintf("{1,6,11}-ROM:\n");

print_mean_time(rom_three.base,"Static Data",exclude_data)
print_mean_time(rom_three.validation_data,"Validation Data",exclude_data)
print_mean_time(rom_three.validation_data-rom_three.base,"Data diff",exclude_data)
print_mean_time(rom_three.orbits,"Orbits",exclude_data)
print_mean_time(rom_three.orbit_validation,"Orbit validation",exclude_data)

rom_three.total = rom_three.validation_data + sum(rom_three.orbits,1) + sum(rom_three.orbit_validation,1);
print_mean_time(rom_three.total,"Total",exclude_data)
fprintf("---\n\n");
fprintf("{1,5,6,11}-ROM:\n");

print_mean_time(rom_four.base,"Static Data",exclude_data)
print_mean_time(rom_four.validation,"Validation Data",exclude_data)
print_mean_time(rom_four.validation-rom_four.base,"Data diff",exclude_data)
print_mean_time(rom_four.orbits,"Orbits",exclude_data)
print_mean_time(rom_four.orbit_validation,"Orbit validation",exclude_data)

rom_four.total = rom_four.validation + sum(rom_four.orbits,1) + sum(rom_four.orbit_validation,1);
print_mean_time(rom_four.total,"Total",exclude_data)
fprintf("---\n\n");

%-----------------------------------
function Static_Data = three_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts)
Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts,"calibration_opts",Calibration_Opts);
Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;
end

function Static_Data = four_mode_rom(system_name,energy_limit,initial_modes,Static_Opts,Verification_Opts,Calibration_Opts)
Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts,"calibration_opts",Calibration_Opts);
Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;
end

%---
%---
%----

function Dyn_Data = one_mode_rom_orbits(system_name)
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 66539;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% 
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e0;
Continuation_Opts.max_inc = 1e0;
Continuation_Opts.min_inc = 1e0;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
% -----------------------------------------%
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
end

function Dyn_Data = one_mode_rom_validation(Dyn_Data)
compare_validation(Dyn_Data,"validation error",1,1:6);
end


function Static_Data = two_mode_rom(Static_Data,added_modes)
Static_Data = Static_Data.update_model(added_modes);
Static_Data = Static_Data.create_dataset;
Static_Data.save_data;
end

function Dyn_Data = two_mode_rom_orbits(system_name,state)
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 66539;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% 
% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e0;
Continuation_Opts.max_inc = 1e0;
Continuation_Opts.min_inc = 1e0;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
% -----------------------------------------%
switch state
    case 1
        Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    case 2
        potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2.69e6,[1e-7,7.5e-8]);
        Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
    case 3
        Continuation_Opts.initial_inc = 5e-2;
        Continuation_Opts.max_inc = 5e-2;
        Continuation_Opts.min_inc = 5e-2;

        Dyn_Data = Dyn_Data.add_orbits(2,[5,8],"opts",Continuation_Opts);
    case 4
        Continuation_Opts.initial_inc = 1e-1;
        Continuation_Opts.max_inc = 1e-1;
        Continuation_Opts.min_inc = 1e-2;
        Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
        Dyn_Data = Dyn_Data.restart_point(2,3,"po","opts",Continuation_Opts);
end
end

function Dyn_Data = two_mode_rom_validation(Dyn_Data,state)
switch state
    case 1
        compare_validation(Dyn_Data,"validation error",[1,2],"all");
    case 2
        compare_validation(Dyn_Data,"validation error",3,[5,11,13]);
    case 3
        Dyn_Data = Dyn_Data.validate_solution(4,[5,11,13]);
        Dyn_Data = Dyn_Data.validate_solution(5,[5,11,13]);
end
end