clear
close all
% start_memory_profiler()
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 0.8;
initial_modes = [1,6];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 1;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 20;
%------------------------------------------%


Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model);
Static_Data.save_data;
% stop_memory_profiler()
% memory_change = get_memory_change;