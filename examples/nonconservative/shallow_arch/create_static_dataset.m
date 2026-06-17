clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "shallow_arch";
%energy_limit = 0.01;
energy_limit = 0.5;
initial_modes = [1,4];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs =  1; %be careful!
% Static_Opts.num_loadcases = 100;
%------------------------------------------%
Verification_Opts.maximum_iterations = 3;


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,"verification_opts",Verification_Opts);
Static_Data.save_data;

% 
% x_bc = Static_Data.get_dataset_values("physical_displacement");
% num_nodes = size(x_bc,1)/3;
% y_index = (0:(num_nodes-1))*6+2;
% 
% y_disp = x_bc(y_index);
% [a,b] = max(y_disp*1000)

%---------------------------------------
