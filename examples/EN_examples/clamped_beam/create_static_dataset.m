clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "clamped_beam";
energy_limit = 0.01;
initial_modes = [1,3];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.output_format = "text";
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs = 1; %be careful!
Static_Opts.num_loadcases = 10;
Static_Opts.maximum_loadcases = 20;
Static_Opts.perturbation_scale_factor = 1;
%------------------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%


%--------- Static Verification Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
Verification_Opts.maximum_iterations = 1;
Verification_Opts.maximum_interpolation_error = [1e-3,1e-3];
Verification_Opts.num_added_points = 1;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;

% Static_Data = load_static_data("clamped_beam_136");
% Static_Opts.max_parallel_jobs = 4; %be careful!
% Static_Opts.additional_data = "stiffness";
% Static_Data = Static_Data.add_additional_data(Static_Opts);
% Static_Data.save_data;