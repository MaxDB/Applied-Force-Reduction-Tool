function create_static_data(initial_modes,additional_data)
%--------- Software Settings ---------%
set_logging_level(1)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 0.8;
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 1;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = additional_data;
Static_Opts.num_validation_modes = 20;
Static_Opts.max_parallel_jobs = 4; %be careful!
%------------------------------------------%


Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model);
Static_Data.save_data;
end