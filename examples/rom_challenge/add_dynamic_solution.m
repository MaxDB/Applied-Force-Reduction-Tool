clear
close all
set_visualisation_level(2)
set_logging_level(2)

system_name = "exhaust_1567";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.type = "physical displacement";
Additional_Output.dof = 1563;
Additional_Output.special_points = [0.25,0.5,0.75,1,1.5,2,3]*1.5e-3;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e1;
Continuation_Opts.max_inc = 1e1;
Continuation_Opts.min_inc = 1e1;
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 2e-1;
% Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.energy_limit_multiplier = 1;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% compare_validation(Dyn_Data,"energy",1,"all")
Dyn_Data = Dyn_Data.validate_solution(1,"all");
% Dyn_Data = Dyn_Data.get_fe_output("periodicity",1,"X");
% Dyn_Data = Dyn_Data.get_max_disp_stress(1,46);

% % 
% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 2e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 500;
% Continuation_Opts.backward_steps = 500;
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.restart_point(2,8,"po","opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.remove_solution(2);
% %-----------------------------------------%


Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 0.67;
Damping_Data.stiffness_factor = 2.55e-7;


Force_Data.type = "point force";
Force_Data.dof = 1563;
Force_Data.continuation_variable = "frequency";
Force_Data.amplitude = 0.17;
% 
%--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e0;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e-1;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% %-----------------------------------------%
% % 

% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% Force_Data.frequency = 1570;
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% % 
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 3e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 1000;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.parameter_range = [0.004,0.0068];
Continuation_Opts.energy_limit_multiplier = 1;
% %-----------------------------------------%
Force_Data.frequency = 942;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);

%--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e0;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e-1;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.restart_point(3,1114,"force","opts",Continuation_Opts);

% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 3e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 500;
% Continuation_Opts.backward_steps = 500;
% Continuation_Opts.parameter_range = [0.004,0.0068];
% Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.restart_point(29,2,"IC");

Dyn_Data = Dyn_Data.get_fe_output("forced_response",3,[1047]);

% 
% 
% Force_Data.continuation_variable = "frequency";
% Force_Data.amplitude = 42;
% bb_sol = [2,73];
% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 5e-2;
% Continuation_Opts.max_inc = 3e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.collation_degree = 6;
% Continuation_Opts.parameter_range = [0.004,0.0068];
% %-----------------------------------------%
% % Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"backbone orbit",bb_sol);
% 
% %-----------------------------------------%
%--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e0;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e-1;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% %-----------------------------------------%
% Force_Data.force_points = [0.1,0.17,0.3];
% Dyn_Data = Dyn_Data.restart_point(3,114,"amplitude");