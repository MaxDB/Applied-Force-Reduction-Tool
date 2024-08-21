clear
close all
set_visualisation_level(3)

system_name = "clamped_beam_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 2e-2;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 8;
Continuation_Opts.energy_limit_multiplier = 1;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% compare_validation(Dyn_Data,1,2:18)
% Dyn_Data = Dyn_Data.validate_solution(1,1:18);




Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 10;
Damping_Data.stiffness_factor = 0;


Force_Data.type = "modal";
Force_Data.mode_number = 1;
Force_Data.continuation_variable = "frequency";
Force_Data.amplitude = 2;
% 

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 3e1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 10;
Continuation_Opts.backward_steps = 500;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.parameter_range = [0.009,0.0209];
%-----------------------------------------%
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);


Dyn_Data = Dyn_Data.get_fe_output("forced_response",2,75);