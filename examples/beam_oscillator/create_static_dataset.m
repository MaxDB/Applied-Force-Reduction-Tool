clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "beam_oscillator";
energy_limit = 1e4;
initial_modes = [1];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.static_solver = "matlab";
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 3;
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_loadcases = 20;
Static_Opts.maximum_loadcases = 25;
%------------------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%

%--------- Static Validation Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
% Validation_Opts.minimum_degree = 3;
Verification_Opts.maximum_iterations = 5;
Verification_Opts.maximum_interpolation_error = [1e-4,1e-4];
Verification_Opts.num_added_points = 1;
Verification_Opts.max_added_points = 800;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);
Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;


