clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mass_spring_axial";

energy_limit = 5e3;
initial_modes = [1];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.max_parallel_jobs = 1; %be careful!
%------------------------------------------%


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model);
Static_Data.save_data;


%----------------------------------------
External_Force.type = "point";
External_Force.dof = 10;

% External_Force.type = "point";
% External_Force.dof = 182;
% External_Force.max_amplitude = 10;

Nc_Data = Nonconservative_Data(Model,External_Force);
Nc_Static_Data = Static_Data.extend_stress_manifold(Nc_Data);

Nc_Static_Data.save_data;
static_dataset_verificiation_plot(Nc_Static_Data)