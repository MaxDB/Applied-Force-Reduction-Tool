clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "clamped_beam";
%energy_limit = 0.01;
energy_limit = 0.007;
initial_modes = [1,3];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "stiffness";
Static_Opts.num_validation_modes = 10;
Static_Opts.max_parallel_jobs =  4; %be careful!
%------------------------------------------%
Verification_Opts.maximum_iterations = 3;


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,"verification_opts",Verification_Opts);
Static_Data.save_data;
%---------------------------------------
if isequal(initial_modes,0)
    return
end
% 
External_Force.type = "point";
External_Force.dof = 248;
External_Force.max_amplitude =  5;

% External_Force.type = "uniform";
% External_Force.direction = 2;
% External_Force.max_amplitude = 10;

% External_Force.type = "point";
% External_Force.dof = 182;
% External_Force.max_amplitude = 10;



% Static_Data = load_static_data("clamped_beam_13");
Nc_Data = Nonconservative_Data(Static_Data.Model,External_Force);
Nc_Static_Data = Static_Data.extend_stress_manifold(Nc_Data);

Nc_Static_Data.save_data;
static_dataset_verification_plot(Nc_Static_Data)




