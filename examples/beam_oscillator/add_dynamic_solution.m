clear
% close all
system_name = "beam_oscillator_1";

Dyn_Data = initalise_dynamic_data(system_name);


% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;

Continuation_Opts.forward_steps = 1500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%
% 
Dyn_Data = Dyn_Data.add_backbone(1,"type","fom","opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.add_backbone(1);
% plot_backbone(system_name,"amplitude","last");
% % compare_validation(Dyn_Data,1,2:18)
% Dyn_Data = Dyn_Data.validate_solution(1,2);
% %plot_h_predicition(system_name,"amplitude",1);
% Dyn_Data = Dyn_Data.add_full_order_backbone(1);
% % Dyn_Data = get_periodicity_error(Dyn_Data,2,1:587)


% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-2;
% Continuation_Opts.max_inc = 1e-2;
% Continuation_Opts.min_inc = 1e-3;
% Continuation_Opts.forward_steps = 300;
% Continuation_Opts.backward_steps = 0;
% % Continuation_Opts.parameter_range = [-10,10];
% %-----------------------------------------%D
% Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);
%  Dyn_Data = Dyn_Data.restart_point(1,199,"po");
% plot_backbone(system_name,"amplitude","last");



% Damping_Data.damping_type = "rayleigh";
% Damping_Data.mass_factor = 0.95;
% Damping_Data.stiffness_factor = 1.58e-6;
% 
% 
% Force_Data.type = "modal";
% Force_Data.mode_number = 1;
% Force_Data.continuation_variable = "amplitude";
% Force_Data.force_points = [10];
% 
% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e0;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e-1;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);
% %-----------------------------------------%
% 
% Force_Data.frequency = 942;
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data);
% Force_Data.frequency = 1570;
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data);
% 
% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 1e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 10;
% Continuation_Opts.backward_steps = 10;
% Continuation_Opts.parameter_range = [0.004,0.0068];
% Dyn_Data = Dyn_Data.update_continuation_opts(Continuation_Opts);
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.restart_point(6,"all","IC");
% Dyn_Data = Dyn_Data.restart_point(3,"all","IC");

% plot_backbone(Dyn_Data,"energy","last");
% plot_backbone(Dyn_Data,"amplitude","last");
% plot_h_predicition(Dyn_Data,"amplitude","last");
% plot_h_predicition(Dyn_Data,"energy","last");

