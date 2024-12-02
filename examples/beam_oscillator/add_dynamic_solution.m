clear
% close all
system_name = "beam_oscillator_1";

Dyn_Data = initalise_dynamic_data(system_name);


% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e1;
Continuation_Opts.max_inc = 1e1;
Continuation_Opts.min_inc = 1e-2;

Continuation_Opts.forward_steps = 1500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%
Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);
Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

Dyn_Data = Dyn_Data.validate_solution(2,2);

