clear
close all
%--------- System Settings ---------%
system_name = "exhaust_157";
added_modes = 6;
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
%------------------------------------------%

load("data\" + system_name + "\" + "Static_Data.mat","Static_Data");
Static_Data = Static_Data.update_model(added_modes,Calibration_Opts,Static_Opts);


%--------- Static Solver Settings ---------%
Static_Opts.num_loadcases = 10;
Static_Opts.maximum_loadcases = 10;
Static_Data.Model = Static_Data.Model.update_static_opts(Static_Opts);
%------------------------------------------%

%--------- Static Validation Settings ---------%
Validation_Opts.validation_algorithm = "sep_points";
% Validation_Opts.minimum_degree = 3;
Validation_Opts.minimum_coupling_rating = 1e-2;
Validation_Opts.maximum_iterations = 3;
Validation_Opts.maximum_interpolation_error = [1e-3,1e-1];
Validation_Opts.maximum_fitting_error = 1e-3;
Validation_Opts.num_added_points = 1;
Validation_Opts.max_added_points = 800;
Static_Data = Static_Data.update_validation_opts(Validation_Opts);
%----------------------------------------------%

Static_Data = Static_Data.create_dataset;
Static_Data.save_data;