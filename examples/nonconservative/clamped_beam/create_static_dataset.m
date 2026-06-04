clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "clamped_beam";
energy_limit = 0.01;
initial_modes = [1];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs =  4; %be careful!
% Static_Opts.num_loadcases = 100;
%------------------------------------------%
Verification_Opts.maximum_iterations = 3;


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,"verification_opts",Verification_Opts);
Static_Data.save_data;



%---------------------------------------

% 
External_Force.type = "point";
External_Force.dof = 236;
External_Force.max_amplitude = 1;

% External_Force.type = "uniform";
% External_Force.direction = 2;
% External_Force.max_amplitude = 10;

% External_Force.type = "point";
% External_Force.dof = 182;
% External_Force.max_amplitude = 10;

set_visualisation_level(3)

Nc_Data = Nonconservative_Data(Model,External_Force);
Nc_Static_Data = Static_Data.extend_stress_manifold(Nc_Data);

Nc_Static_Data.save_data;
static_dataset_verification_plot(Nc_Static_Data)


%------
Rom = Reduced_System(Static_Data);
Rom_Nc = Reduced_System(Nc_Static_Data);


