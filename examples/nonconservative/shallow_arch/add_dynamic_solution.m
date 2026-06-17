clear
close all
set_visualisation_level(1)

system_name = "shallow_arch_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [0.32,0.00192,0.016];
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 3e-1;
Continuation_Opts.max_inc = 3e-1;
Continuation_Opts.min_inc = 1e-2;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
%-----------------------------------------%

% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

%-------------------
%Forced response
%-------------------
Continuation_Opts.collocation_degree = 10;
Continuation_Opts.min_discretisation_num = 40;
Continuation_Opts.initial_discretisation_num = 40;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.frequency_range = [1e6,1.04e6];

%-----------------------------------------%
Model = Dyn_Data.Dynamic_Model.Model;
Force_Data.type = "shape";
Force_Data.frequency = 1.015e6;
Force_Data.amplitude = 1;
Force_Data.continuation_variable = "frequency";

force_unit_factor = 1e3;
kappa_1 = 0.03*force_unit_factor; %(mm/s^2)


%---
Damping_Data.damping_type = "nonlinear_rayleigh";
damping_coeffs = get_shallow_arch_damping(Model);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);
%---

kappa_2 = 10*force_unit_factor;
Force_Data.shape = get_shallow_arch_force(Model,kappa_1,kappa_2);
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
return

kappa_2 = -10e3;
Force_Data.shape = get_shallow_arch_force(Model,kappa_1,kappa_2);
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);



function force_shape = get_shallow_arch_force(Model,kappa_1,kappa_2)
mass = Model.mass;
stiffness = Model.stiffness;

[evec,~] = eigs(stiffness,mass,4,"smallestabs");
bending_13 = evec(:,[1,4]);
force_shape = mass*bending_13*[kappa_1;kappa_2];
end

function damping_coeffs = get_shallow_arch_damping(Model)
mass = Model.mass;
stiffness = Model.stiffness;

[~,eval] = eigs(stiffness,mass,1,"smallestabs");

freq_1 = sqrt(eval);

damping_coeffs(1) = freq_1/500;
damping_coeffs(2) = 0;
end