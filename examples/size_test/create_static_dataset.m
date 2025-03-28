clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 1;
initial_modes = [1];
%-----------------------------------%

%--------- Calibration Settings ---------%
% Calibration_Opts.Static_Opts.num_loadcases = 20;
Calibration_Opts.Static_Opts = struct([]);
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.output_format = "binary";
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.num_added_points = 1;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;