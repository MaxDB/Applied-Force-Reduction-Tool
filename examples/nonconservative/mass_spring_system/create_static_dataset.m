clear
close all
%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(3)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "mass_spring_roller";

% energy_limit = 0.05;
energy_limit = 3e-3;
initial_modes = [1,2];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.additional_data = "stiffness";
Static_Opts.max_parallel_jobs = 1; %be careful!

%------------------------------------------%


Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model);
Static_Data.save_data;
%----------------------------------------

return

External_Force.type = "shape";
External_Force.shape = get_angled_force_shape(1.8,"deg");
External_Force.max_amplitude = [];

Nc_Data = Nonconservative_Data(Model,External_Force);
Nc_Static_Data = Static_Data.extend_stress_manifold(Nc_Data);
Nc_Static_Data.save_data;


Rom = Reduced_System(Nc_Static_Data);

ax = plot_static_data("force",Nc_Static_Data);
Rom.Force_Polynomial.plot_polynomial("axes",ax);
%
ax = plot_static_data("displacement",Nc_Static_Data);
Rom.Physical_Displacement_Polynomial.plot_polynomial("axes",ax);
%
ax = plot_static_data("energy",Nc_Static_Data);
Rom.Potential_Polynomial.plot_polynomial("axes",ax);


function force_shape = get_angled_force_shape(angle,type)
if type == "deg"
    angle = angle*pi/180;
end
force_shape = [cos(angle);sin(angle);0];
force_shape = force_shape/norm(force_shape);
end