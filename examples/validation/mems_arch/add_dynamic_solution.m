clear

close all
set_visualisation_level(1)
set_logging_level(3)

system_name = "mems_arch_16";
Dyn_Data = initalise_dynamic_data(system_name);
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [0,36e-3,10e-3];
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 5e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 40;
% -----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2.69e6,[1.5e-7,1e-7]);
Continuation_Opts.collocation_degree = 10;
Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);
