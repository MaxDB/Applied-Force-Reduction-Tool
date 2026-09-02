clear
% close all
set_visualisation_level(12)

system_name = "mass_spring_roller_1";
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
Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

%-----------------------------------
%--- angled force
%-----------------------------------
Continuation_Opts.initial_inc = 2e-2;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 5e-3;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.frequency_range = [1.3,2];
%-----------------------------------------%

Force_Data.type = "shape";
Force_Data.shape = get_force_shape(1.8,"deg");
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 1.9;


Damping_Data.damping_type = "rayleigh";
Damping_Data.stiffness_factor = 0;

%---
Force_Data.amplitude = 0.07;
target_damping = 0.01;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,0]);
Damping_Data.mass_factor = damping_coeffs(1);

% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%---
Force_Data.amplitude = 0.7;
target_damping = 0.1;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,0]);
Damping_Data.mass_factor = damping_coeffs(1);

% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%---
Force_Data.amplitude = 3;
target_damping = 1;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,0]);
Damping_Data.mass_factor = damping_coeffs(1);
% 
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%---




function force_shape = get_force_shape(angle,type)
if type == "deg"
    angle = angle*pi/180;
end
force_shape = [cos(angle);sin(angle);0];
force_shape = force_shape/norm(force_shape);
end


%------
%validation
Model = Dyn_Data.Dynamic_Model.Model;

External_Force.type = "shape";
External_Force.shape = get_force_shape(1.8,"deg");
External_Force.max_amplitude = [];


Dyn_Data = Dyn_Data.add_nc_validation_shape(Force_Data);


Dyn_Data = Dyn_Data.validate_solution(1,1001);