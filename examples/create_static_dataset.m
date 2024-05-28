clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust";
% energy_limit = 1.5;
% energy_limit = 0.005;
energy_limit = 6;
initial_modes = [1];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
% Calibration_Opts.energy_overfit = 1.1;
% Calibration_Opts.force_overcalibration = 1.1;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.static_solver = "abaqus";
Static_Opts.additional_data = "perturbation";
Static_Opts.num_validation_modes = 7;
Static_Opts.num_fe_cpus = 1; %not that useful
Static_Opts.max_parallel_jobs = 1; %be careful!
% Static_Opts.num_loadcases = 20;
% Static_Opts.maximum_step_increments = 100;
% Static_Opts.initial_time_increment = 1;
% Static_Opts.total_step_time = 1;
% Static_Opts.minimum_time_increment = 0;
% Static_Opts.maximum_time_increment = 1;
%------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

%--------- Static Solver Settings ---------%
Static_Opts.num_loadcases = 10;
Static_Opts.maximum_loadcases = 10;
Model = Model.update_static_opts(Static_Opts);
%------------------------------------------%

%--------- Static Validation Settings ---------%
Validation_Opts.validation_algorithm = "sep_points";
% Validation_Opts.minimum_degree = 3;
Validation_Opts.minimum_coupling_rating = 1e-2;
Validation_Opts.maximum_iterations =3;
Validation_Opts.maximum_interpolation_error = [1e-2,1e-1];
Validation_Opts.maximum_fitting_error = 1e-4;
Validation_Opts.num_added_points = 1;
%----------------------------------------------%
Static_Data = Static_Dataset(Model,Validation_Opts);
Static_Data.save_data;

