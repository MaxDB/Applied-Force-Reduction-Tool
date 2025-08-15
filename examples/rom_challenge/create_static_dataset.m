clear
close all
%--------- Software Settings ---------%
set_logging_level(4)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust";
energy_limit = 1.8;
initial_modes = [1,7];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.num_validation_modes = 18;
Static_Opts.max_parallel_jobs = 4; %be careful!
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
Verification_Opts.maximum_iterations =3;
Verification_Opts.num_added_points =1;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;