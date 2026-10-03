clear
close all
%--------- Software Settings ---------%
set_logging_level(4)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust_1";
added_modes = 7;
%-----------------------------------%

%--------- Static Validation Settings ---------%
Verification_Opts.maximum_iterations = 8;
Verification_Opts.num_added_points = 4;


% Validation_Opts.maximum_iterations = 3;
% Validation_Opts.max_added_points = 150;

% Validation_Opts.maximum_iterations = 2;
% Validation_Opts.max_added_points = 800;
%------------------------------------------%

Static_Data = load_static_data(system_name);
Static_Data = Static_Data.update_verification_opts(Verification_Opts);
Static_Data = Static_Data.update_model(added_modes);
Static_Data = Static_Data.create_dataset;
Static_Data.save_data;