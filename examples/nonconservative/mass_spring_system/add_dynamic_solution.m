clear
% close all
set_visualisation_level(1)

system_name = "mass_spring_roller_0";
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 2;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;

Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 8;


%-----------------------------------------%
Continuation_Opts.energy_limit_multiplier = 1.2;
% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

%-----------------------------------
%--- Vertical force
%-----------------------------------
Continuation_Opts.initial_inc = 2e-2;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 5e-3;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.frequency_range = [1.3,3.2];
%-----------------------------------------%

Force_Data.type = "shape";
Force_Data.shape = [0;1;1];
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 1.35;
% Damping_Data.damping_type = "nonlinear_rayleigh";
Damping_Data.damping_type = "rayleigh";

% Force_Data.amplitude = 0.001;
Force_Data.amplitude = 0.008;
target_damping = 0.01;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%---

Force_Data.amplitude = 0.002;
target_damping = 0.2;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);

Force_Data.amplitude = 0.09;
target_damping = 1;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);