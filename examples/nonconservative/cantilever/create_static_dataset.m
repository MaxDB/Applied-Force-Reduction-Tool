clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "cantilever_f";
energy_limit = 0.76;
initial_modes = [2001];
%-----------------------------------%

% Calibration_Opts.calibration_scale_factor = 2;
Static_Opts.num_loadcases = 20;

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.max_parallel_jobs = 1; %be careful!
% Static_Opts.follower_validation_modes = 1;
%------------------------------------------%

% Verification_Opts.maximum_iterations = -1;
Verification_Opts.maximum_iterations = 3;

Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,"verification_opts",Verification_Opts);

%---
Static_Data.verified_degree = [9,9];
static_dataset_verification_plot(Static_Data)

%---

Static_Data.save_data;