% clear
% close all
%--------- System Settings ---------%
system_name = "h_oscillator_1";
added_modes = 2;
%-----------------------------------%

load("data\" + system_name + "\" + "Static_Data.mat","Static_Data");
% Static_Data = Static_Data.update_model(added_modes,Static_Opts,Calibration_Opts);
Static_Data = Static_Data.update_model(added_modes);

Static_Data = Static_Data.create_dataset;
Static_Data.save_data;