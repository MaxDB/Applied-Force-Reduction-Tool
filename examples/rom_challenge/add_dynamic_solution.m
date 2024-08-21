clear
close all
set_visualisation_level(3)

system_name = "exhaust_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.type = "physical displacement";
Additional_Output.dof = 1563;
Additional_Output.special_points = [0.25,0.5,0.75,1,1.5,2,3]*1.5e-3;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e1;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e1;
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.energy_limit_multiplier = 1;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% compare_validation(Dyn_Data,1,2:18)
% Dyn_Data = Dyn_Data.validate_solution(1,1:18);
% Dyn_Data = Dyn_Data.get_fe_output("periodicity",1,"X");
% Dyn_Data = Dyn_Data.get_periodicity_error(1,"all");
% Dyn_Data = Dyn_Data.get_periodicity_error(1,"X");
% Dyn_Data = Dyn_Data.get_max_disp_stress(1,46);


% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 2e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 500;
% Continuation_Opts.backward_steps = 500;
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.restart_point(2,8,"po","opts",Continuation_Opts);
% %-----------------------------------------%


Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 0.95;
Damping_Data.stiffness_factor = 1.58e-6;


Force_Data.type = "modal";
Force_Data.mode_number = 1;
Force_Data.continuation_variable = "amplitude";
Force_Data.force_points = [10,50,100,150,200,250,300,350,400];
% 
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e0;
Continuation_Opts.max_inc = 1e1;
Continuation_Opts.min_inc = 1e-1;
Continuation_Opts.forward_steps = 100;
Continuation_Opts.backward_steps = 0;
%-----------------------------------------%

Force_Data.frequency = 942;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
Force_Data.frequency = 1570;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% % 
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 3e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 100;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.parameter_range = [0.004,0.0068];
%-----------------------------------------%
Dyn_Data = Dyn_Data.restart_point(2,6,"IC","opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.restart_point(4,2:9,"IC");
% 
% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e0;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e-1;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.restart_point(21,132,"force");
% 
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
