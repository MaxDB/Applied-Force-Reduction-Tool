clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust";
energy_limit = 1.8;
initial_modes = [1,7];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.static_solver = "abaqus";
Static_Opts.additional_data = "perturbation";
Static_Opts.num_validation_modes = 18;
Static_Opts.max_parallel_jobs = 1; %be careful!
Static_Opts.num_loadcases = 8;
Static_Opts.maximum_loadcases = 15;
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
Verification_Opts.maximum_iterations = 3;
Verification_Opts.maximum_interpolation_error = [1e-3,1e-2];
Verification_Opts.num_added_points = 3;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;