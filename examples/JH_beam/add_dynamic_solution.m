clear
% close all
set_visualisation_level(1)

system_name = "JH_beam_2d_135";
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
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 1e-2;
Continuation_Opts.forward_steps = 000;
Continuation_Opts.backward_steps = 250;
Continuation_Opts.initial_discretisation_num = 200;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(1,1,"po","opts",Continuation_Opts);





%-------------------------------------------
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 1e-2;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 250;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(1,11,"po","opts",Continuation_Opts);

%-------------------------------------------
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 1e-2;
Continuation_Opts.forward_steps = 300;
Continuation_Opts.backward_steps = 300;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(1,17,"po","opts",Continuation_Opts);

%-------------------------------------------
Continuation_Opts.initial_inc = 5e-3;
Continuation_Opts.max_inc = 5e-3;
Continuation_Opts.forward_steps = 00;
Continuation_Opts.backward_steps = 500;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(3,218,"po","opts",Continuation_Opts);


Continuation_Opts.initial_inc = 5e-3;
Continuation_Opts.max_inc = 5e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 500;
%-----------------------------------------
Dyn_Data = Dyn_Data.restart_point(3,199,"po","opts",Continuation_Opts);

