clear
% close all
set_visualisation_level(1)

system_name = "mass_spring_roller_12";
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 2;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 2e-1;
% Continuation_Opts.max_inc = 2e-2;

Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 2000;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 100;
Continuation_Opts.max_discretisation_num = 250;

Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 8;
Continuation_Opts.energy_limit_multiplier = 1.2;
%-----------------------------------------%
 % Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);
Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

Dyn_Data = Dyn_Data.validate_solution(1,2);

%c
% compare_validation(Dyn_Data,"validation error",1,"all")


% compare_validation(Dyn_Data,"force amplitude","last","all")
% compare_solutions("energy",Dyn_Data,1)

% 
% plot_h_predicition(Dyn_Data,"mean error",1);