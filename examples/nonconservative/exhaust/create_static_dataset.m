clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "exhaust";
energy_limit = 1;
initial_modes = [1,7];
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 1.5;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.num_validation_modes = 18;
Static_Opts.max_parallel_jobs = 1; %be careful!
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.maximum_iterations = 3;
%----------------------------------------------%

Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);

Static_Data = Static_Dataset(Model,"verification_opts",Verification_Opts);
Static_Data.save_data;
%----------------------------------------------%
if isequal(initial_modes,0)
    return
end

External_Force.type = "uniform";
External_Force.direction =3;
External_Force.max_amplitude = 10;


set_visualisation_level(3)

Nc_Data = Nonconservative_Data(Static_Data.Model,External_Force);
Nc_Static_Data = Static_Data.extend_stress_manifold(Nc_Data);

Nc_Static_Data.save_data;
static_dataset_verification_plot(Nc_Static_Data)
