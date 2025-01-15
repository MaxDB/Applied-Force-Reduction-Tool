clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 2;
initial_modes = [1];
%-----------------------------------%

% Add mesh data file so dimension can be extracted when needed
%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.static_solver = "abaqus";
Static_Opts.additional_data = "none";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_loadcases = 10;
Static_Opts.maximum_loadcases = 20;
%------------------------------------------%

%--------- Static Validation Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
Verification_Opts.maximum_iterations = 3;
Verification_Opts.maximum_interpolation_error = [1e-3,1e-2];
Verification_Opts.num_added_points = 1;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;