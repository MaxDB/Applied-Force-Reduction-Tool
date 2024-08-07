clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust_17";
added_modes = 5;
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.num_loadcases = 6;
Static_Opts.maximum_loadcases = 9;
%------------------------------------------%

%--------- Static Validation Settings ---------%
% Validation_Opts.maximum_iterations = 3;
% Validation_Opts.max_added_points = 50;
% 
Validation_Opts.maximum_iterations = 3;
Validation_Opts.max_added_points = 150;

% Validation_Opts.maximum_iterations = 2;
% Validation_Opts.max_added_points = 800;
%------------------------------------------%

load("data\" + system_name + "\static_data\" + "Static_Data.mat","Static_Data");
Static_Data = Static_Data.update_validation_opts(Validation_Opts);
Static_Data = Static_Data.update_model(added_modes,Static_Opts);
Static_Data = Static_Data.create_dataset;
Static_Data.save_data;