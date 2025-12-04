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
initial_mode = 1;
added_mode = 6;
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "stiffness";
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_validation_modes = 20;
%------------------------------------------%
Calibration_Opts.calibration_scale_factor = 1;
%------------------------------------------%




rom_one_creation = zeros(1,num_iterations);
rom_one_orbits = zeros(1,num_iterations);
rom_one_orbit_validation = zeros(1,num_iterations);

rom_two_creation = zeros(1,num_iterations);
rom_two_orbits = zeros(2,num_iterations);
rom_two_orbit_validation = zeros(2,num_iterations);


num_dof = create_mesh(0.00307);

create_parallel_pool(Static_Opts.max_parallel_jobs);
for iCount = 1:num_iterations
    close all
    delete_static_data(system_name+"_"+initial_mode);
    delete_cache(system_name,"force",energy_limit)
    delete_cache(system_name,"matrices")

    validation_time_start = tic;
    Static_Data = one_mode_rom(system_name,energy_limit,initial_mode,Static_Opts,Calibration_Opts);
    rom_one_creation(iCount) = toc(validation_time_start);

    orbit_time_start = tic;
    Dyn_Data = one_mode_rom_orbits(system_name+"_"+initial_mode);
    rom_one_orbits(1,iCount) = toc(orbit_time_start);

    orbit_validation_start = tic;
    Dyn_Data = one_mode_rom_validation(Dyn_Data);
    rom_one_orbit_validation(1,iCount) = toc(orbit_validation_start);


    clear("Dyn_Data")
    close all
    % %----------------------------------------------------------------
    validation_time_start = tic;
    Static_Data = two_mode_rom(Static_Data,added_mode);
    rom_two_creation(1,iCount) = toc(validation_time_start);
    %

    two_mode_name = system_name+"_"+initial_mode + added_mode;
    for iOrbit = 1:2
        orbit_time_start = tic;
        Dyn_Data = two_mode_rom_orbits(two_mode_name,iOrbit);
        rom_two_orbits(iOrbit,iCount) = toc(orbit_time_start);
    end

    %
    for iValidation = 1:2
        orbit_validation_start = tic;
        Dyn_Data = two_mode_rom_validation(Dyn_Data,iValidation);
        rom_two_orbit_validation(iValidation,iCount) = toc(orbit_validation_start);
    end
    %----------------------------------------------------------------
    clear("Dyn_Data")
    clear("Static_Data")

end
exclude_data = [];
fprintf("\n\n1-ROM:\n");

print_mean_time(rom_one_creation,"Static Data",exclude_data)
print_mean_time(rom_one_orbits,"Orbits",exclude_data)
print_mean_time(rom_one_orbit_validation,"Orbit validation",exclude_data)

total_one = rom_one_creation + sum(rom_one_orbits,1) + sum(rom_one_orbit_validation,1);
print_mean_time(total_one,"Total",exclude_data)
fprintf("---\n\n");
fprintf("{1,6}-ROM:\n");

print_mean_time(rom_two_creation,"Static Data",exclude_data)
print_mean_time(rom_two_orbits(1,:),"Orbits 1",exclude_data)
print_mean_time(rom_two_orbits(2,:),"Orbits 2",exclude_data)
print_mean_time(rom_two_orbit_validation(1,:),"Orbit validation 1",exclude_data)
print_mean_time(rom_two_orbit_validation(2,:),"Orbit validation 2",exclude_data)

total_two = rom_two_creation + sum(rom_two_orbits,1) + sum(rom_two_orbit_validation,1);
print_mean_time(total_two,"Total",exclude_data)
fprintf("---\n\n");
%-----------------------------------
function Static_Data = one_mode_rom(system_name,energy_limit,modes,Static_Opts,Calibration_Opts)
Model = Dynamic_System(system_name,energy_limit,modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model);
Static_Data.save_data;
end

function Dyn_Data = one_mode_rom_orbits(system_name)
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
switch state
    case 1
        Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    case 2
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
        Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
end

end

function Dyn_Data = two_mode_rom_validation(Dyn_Data,state)

switch state
    case 1
        compare_validation(Dyn_Data,"validation error",[1,2],"all")
    case 2
        Dyn_Data.validate_solution(1,[5,11,13]);
        Dyn_Data.validate_solution(2,[5,11,13],"load_data",1);
end
end