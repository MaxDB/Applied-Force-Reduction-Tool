clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "cantilever_f";
energy_limit = 0.5;
initial_modes = [1];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.max_parallel_jobs = 1; %be careful!
Static_Opts.follower_force = 1;
Static_Opts.output_format = "text";
%------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model);
Static_Data.save_data;