clear
% close all
set_visualisation_level(1)
set_logging_level(2)

system_name = "cubic_mass_spring_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 1;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 8;
Continuation_Opts.energy_limit_multiplier = 2;
%-----------------------------------------%
Dyn_Data = Dyn_Data.add_backbone(1,"type","rom","opts",Continuation_Opts);
Continuation_Opts.inertial_compensation = 0;
Dyn_Data = Dyn_Data.add_backbone(1,"type","rom","opts",Continuation_Opts);
%-----
Continuation_Opts.energy_limit_multiplier = 30;
Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);


Continuation_Opts.initial_inc = 1e0;
Continuation_Opts.max_inc = 1e0;
Continuation_Opts.min_inc = 1e0;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 5;
Dyn_Data = Dyn_Data.restart_point(3,33,"po","opts",Continuation_Opts);

%-----
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 500;
Dyn_Data = Dyn_Data.restart_point(4,5,"po","opts",Continuation_Opts);

%---
Dyn_Data = Dyn_Data.remove_solution(4);

%--




%-----------------------------------------%
% Continuation_Opts.energy_limit_multiplier = 100;
% Continuation_Opts.initial_inc = 3e0;
% Continuation_Opts.max_inc = 3e0;
% Continuation_Opts.min_inc = 3e0;
% Continuation_Opts.forward_steps = 10;
% Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);
% %-----------------------------------------%
% 
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 2e-1;
% Continuation_Opts.min_inc = 1e-2;
% 
% Continuation_Opts.forward_steps = 300;
% Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);
% 
% 
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 100;
% Dyn_Data = Dyn_Data.restart_point(1,3,"po","opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.restart_point(1,5,"po","opts",Continuation_Opts);
