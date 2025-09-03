clear
close all
%--------- Software Settings ---------%
set_logging_level(4)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "JH_beam_2d";
energy_limit = 0.015; %0.01 J
initial_modes = [1,3,5];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs =  4; %be careful!
%------------------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts = struct([]);
%----------------------------------------%


%--------- Static Verification Settings ---------%
Verification_Opts.maximum_iterations = 3;
% [1e-3] works for three modes and 5e-3 works for two modes]
%----------------------------------------------%W
tic
Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,"verification_opts",Verification_Opts);
Static_Data.save_data;
toc
