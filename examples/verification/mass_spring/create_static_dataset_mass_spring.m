clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "cubic_mass_spring";
energy_limit = 0.04;
initial_modes = [1,2];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 1;
Calibration_Opts.Static_Opts.num_loadcases = 300;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.solver_algorithm = "riks";
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_loadcases = 200;
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 2;
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
Verification_Opts.maximum_iterations = 0;
Verification_Opts.maximum_interpolation_error = [1e-3,1e-2];
Verification_Opts.num_added_points = 1;
%----------------------------------------------%
Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;