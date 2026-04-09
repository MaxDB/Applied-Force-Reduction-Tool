clear
close all
set_visualisation_level(2)
set_logging_level(2)

system_name = "cantilever_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
% Additional_Output.type = "physical displacement";
% Additional_Output.dof = 722;
% Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 6;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
%-----------------------------------------%
