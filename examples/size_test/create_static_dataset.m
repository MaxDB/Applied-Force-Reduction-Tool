clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mems_arch";
energy_limit = 0.8;
initial_modes = [1,6];
%-----------------------------------%

%--------- Calibration Settings ---------%
% Calibration_Opts.Static_Opts.num_loadcases = 20;
Calibration_Opts.Static_Opts = struct([]);
Calibration_Opts.calibration_scale_factor = 1.5;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 20;
Static_Opts.output_format = "binary";
% Static_Opts.num_loadcases = 3;
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.num_added_points = 1;
Verification_Opts.maximum_iterations = 1;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;

% Static_Data = load_static_data("mems_arch_1611");
% Static_Opts.additional_data = "perturbation";
% Static_Opts.num_validation_modes = 14;
% Static_Opts.output_format = "text";
% Static_Opts.max_parallel_jobs = 1; %be careful!
% Static_Data = Static_Data.add_additional_data(Static_Opts);
% Static_Data.save_data;