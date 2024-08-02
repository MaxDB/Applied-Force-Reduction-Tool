clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust";
energy_limit = 1.8;
initial_modes = [1];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.static_solver = "abaqus";
Static_Opts.additional_data = "perturbation";
Static_Opts.num_validation_modes = 18;
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_loadcases = 8;
Static_Opts.maximum_loadcases = 15;
%------------------------------------------%

%--------- Static Validation Settings ---------%
Validation_Opts.validation_algorithm = "sep_points_new";
Validation_Opts.minimum_coupling_rating = 1e-2;
Validation_Opts.maximum_iterations = 0;
Validation_Opts.maximum_interpolation_error = [1e-3,1e-3];
Validation_Opts.maximum_fitting_error = 1e-3;
Validation_Opts.num_added_points = 1;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

Static_Data = Static_Dataset(Model,Validation_Opts);
Static_Data.save_data;