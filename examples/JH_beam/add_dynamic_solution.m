clear
% close all
set_visualisation_level(1)

system_name = "JH_beam_2d_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 5e-2;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 200;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% compare_validation(Dyn_Data,"validation error",1,"all")


%-------------------------------------------
Continuation_Opts.initial_inc = 1e-3;
Continuation_Opts.max_inc = 5e-3;
Continuation_Opts.forward_steps = 250;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 200;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(2,1,"po","opts",Continuation_Opts);





%-------------------------------------------
Continuation_Opts.initial_inc = 5e-3;
Continuation_Opts.max_inc = 5e-3;
Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 1000;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(2,255,"po","opts",Continuation_Opts);


%-------------------------------------------
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 1e-2;
Continuation_Opts.forward_steps = 250;
Continuation_Opts.backward_steps = 250;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(1,137,"po","opts",Continuation_Opts);


%-------------------------------------------
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 2e-2;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 128;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(2,252,"po","opts",Continuation_Opts);