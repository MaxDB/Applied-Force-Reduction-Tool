clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mass_spring_roller";

energy_limit = 0.06;
initial_modes = [1];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "none";
Static_Opts.max_parallel_jobs = 4; %be careful!
%------------------------------------------%


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model);
Static_Data.save_data;


%----------------------------------------
% Damping_Data.damping_type = "rayleigh";
% Damping_Data.mass_factor = 0.67;
% Damping_Data.stiffness_factor = 2.55e-7;


% Force_Data.type = "point force";
% Force_Data.dof = 1563;
% Force_Data.continuation_variable = "frequency";
% Force_Data.amplitude = 0.17;
% Force_Data.frequency = 942;

% Force_Data.type = "modal";
% Force_Data.mode_number = 1;
% Force_Data.continuation_variable = "amplitude";
% Force_Data.force_points = [1,10,50];



Damping.type = "rayleigh";
Damping.mass_factor = 0.1;
Damping.stiffness_factor = 1e-7;

External_Force.type = "point";
External_Force.dof = 3;
External_Force.max_amplitude = 1; %"limit" -> calibrate

Nc_Data = Nonconservative_Data(Model,Damping,External_Force);
Nc_Static_Data = Static_Data.extend_stress_manifold(Nc_Data);

Nc_Static_Data.verified_degree = [7,7];


Rom = Reduced_System(Nc_Static_Data);

ax = plot_static_data("force",Nc_Static_Data);
Rom.Force_Polynomial.plot_polynomial("axes",ax);
%
ax = plot_static_data("displacement",Nc_Static_Data);
Rom.Physical_Displacement_Polynomial.plot_polynomial("axes",ax);
%
ax = plot_static_data("energy",Nc_Static_Data);
Rom.Potential_Polynomial.plot_polynomial("axes",ax);

animate_whisker(Rom.Physical_Displacement_Polynomial.subpoly(3),Nc_Data,3,"energy",{Rom.Potential_Polynomial,Model.energy_limit})
animate_whisker(Rom.Force_Polynomial.subpoly(1),Nc_Data,3,"energy",{Rom.Potential_Polynomial,Model.energy_limit})
animate_whisker(Rom.Potential_Polynomial,Nc_Data,3,"energy",{Rom.Potential_Polynomial,Model.energy_limit})