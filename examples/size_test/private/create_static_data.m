function create_static_data(initial_modes,additional_data)
%--------- Software Settings ---------%
set_logging_level(1)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 2;
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.Static_Opts = struct([]);
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = additional_data;
Static_Opts.num_validation_modes = 19;
Static_Opts.max_parallel_jobs = 4; %be careful!
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.num_added_points = 1;
%----------------------------------------------%


Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;
end