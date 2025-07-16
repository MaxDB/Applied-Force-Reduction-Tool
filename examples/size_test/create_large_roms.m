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
initial_modes = [1,6];
added_modes = [11,5];
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

if isempty(gcp('nocreate')) && Static_Opts.max_parallel_jobs > 1
    parpool;
end

Static_Data_16 = load_static_data("mems_arch_16");
Static_Data_16 = Static_Data_16.update_verification_opts(Verification_Opts);


rom_three_name = "mems_arch_1611";
rom_three_base = zeros(1,num_iterations);
rom_three_validation_data = zeros(1,num_iterations);
rom_three_orbits = zeros(1,num_iterations);
rom_three_orbit_validation = zeros(1,num_iterations);

rom_four_base = zeros(1,num_iterations);
rom_four_validation = zeros(1,num_iterations);
rom_four_orbits = zeros(4,num_iterations);
rom_four_orbit_validation = zeros(3,num_iterations);



for iCount = 1:num_iterations
    close all
    delete_static_data(rom_three_name );
    % delete_cache(system_name,"force",energy_limit) %only delete three mode ROM
    
    base_time_start = tic;
    Static_Opts.additional_data = "none";
    Static_Data_Base = three_mode_rom(Static_Data_16,added_modes(1));
    rom_three_base(iCount) = toc(base_time_start);

   
   
end
exclude_data = [];
fprintf("\n\n1-ROM:\n");

print_mean_time(rom_three_base,"Static Data",exclude_data)
print_mean_time(rom_three_validation_data,"Validation Data",exclude_data)
print_mean_time(rom_three_validation_data-rom_three_base,"Data diff",exclude_data)
print_mean_time(rom_three_orbits,"Orbits",exclude_data)
print_mean_time(rom_three_orbit_validation,"Orbit validation",exclude_data)

total_one = rom_three_validation_data + sum(rom_three_orbits,1) + sum(rom_three_orbit_validation,1);
print_mean_time(total_one,"Total",exclude_data)
fprintf("---\n\n");
fprintf("{1,6}-ROM:\n");

print_mean_time(rom_four_base,"Static Data",exclude_data)
print_mean_time(rom_four_validation,"Validation Data",exclude_data)
print_mean_time(rom_four_validation-rom_four_base,"Data diff",exclude_data)
print_mean_time(rom_four_orbits,"Orbits",exclude_data)
print_mean_time(rom_four_orbit_validation,"Orbit validation",exclude_data)

total_two = rom_four_validation + sum(rom_four_orbits,1) + sum(rom_four_orbit_validation,1);
print_mean_time(total_two,"Total",exclude_data)
fprintf("---\n\n");

%-----------------------------------
function Static_Data = three_mode_rom(Static_Data,added_modes)
Static_Data = Static_Data.update_model(added_modes);
Static_Data = Static_Data.create_dataset;
Static_Data.save_data;
end

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