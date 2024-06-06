clear
% close all
system_name = "exhaust_157";
Additional_Output.type = "physical displacement";
Additional_Output.dof = 1563;
Additional_Output.special_points = [0.25,0.5,0.75,1,1.5,2,3]*1.5e-3;

%-------------------------------------------------------------------------%
Dyn_Data = initalise_dynamic_data(system_name,Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.initial_inc = 1e1;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e1;
Continuation_Opts.forward_steps = 1500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.energy_limit_multiplier = 1;
%-----------------------------------------%
Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);

Dyn_Data = Dyn_Data.add_backbone(1);
% compare_validation(Dyn_Data,1,2:18)
% Dyn_Data = Dyn_Data.validate_solution(1,2);
%
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 500;
%-----------------------------------------%D
Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);
 Dyn_Data = Dyn_Data.restart_backbone(1,9,"po");


% Dyn_Data = Dyn_Data.add_full_order_backbone(1);

% Damping_Data.damping_type = "rayleigh";
% Damping_Data.alpha = 4;
% Damping_Data.beta = 1e-6;

% Force_Data.amplitude = [1;0];
% Dyn_Data = Dyn_Data.add_modal_frf(Force_Data,Damping_Data);

% plot_backbone(Dyn_Data,"energy","last");
% plot_backbone(Dyn_Data,"amplitude","last");
% plot_h_predicition(Dyn_Data,"amplitude","last");
% plot_h_predicition(Dyn_Data,"energy","last");

