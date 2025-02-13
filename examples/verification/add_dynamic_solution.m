clear
% close all
set_visualisation_level(1)
set_logging_level(2)

system_name = "ic_demo_12";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "amplitude";
Additional_Output.dof = 362;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
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
Continuation_Opts.inertial_compensation = 1;
Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
Continuation_Opts.inertial_compensation = 0;
Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);
% compare_validation(Dyn_Data,"energy",1,"all")
% Dyn_Data = Dyn_Data.validate_solution(1,"all");
% 
% plot_h_predicition(Dyn_Data,"amplitude",1,"backbone",0);

