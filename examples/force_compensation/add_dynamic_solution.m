clear
close all
set_visualisation_level(2)
set_logging_level(2)

system_name = "fc_cantilever_12";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.type = "physical displacement";
Additional_Output.dof = 722;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
%-----------------------------------------%
Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 4;
Damping_Data.stiffness_factor = 0;


Force_Data.type = "point force";
Force_Data.dof = 722;
Force_Data.continuation_variable = "frequency";
Force_Data.amplitude = 0.9;
Force_Data.frequency = 45;
% 
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
Continuation_Opts.parameter_range = [0.089,0.14];
Continuation_Opts.energy_limit_multiplier = 1;
% %-----------------------------------------%
% 
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% % 
% 
 % Dyn_Data = Dyn_Data.get_fe_output("forced_response",2,82);
