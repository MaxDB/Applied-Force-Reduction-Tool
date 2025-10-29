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
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs =  1; %be careful!
%------------------------------------------%



Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model);
Static_Data.save_data;
