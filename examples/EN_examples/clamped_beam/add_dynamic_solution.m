clear
% close all
set_visualisation_level(1)

system_name = "clamped_beam_13";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 5e-2;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 8;
Continuation_Opts.energy_limit_multiplier = 1;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
compare_validation(Dyn_Data,"validation error",1,3)
% Dyn_Data = Dyn_Data.validate_solution(1,6);
Continuation_Opts.initial_inc = 1e-3;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.max_inc = 1e-4;
Continuation_Opts.parameter_range = [0.0146,100];
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 200;
Dyn_Data = Dyn_Data.restart_point(1,4,"po","opts",Continuation_Opts);
% 
% compare_validation(Dyn_Data,"validation error","last",6)
% Damping_Data.damping_type = "rayleigh";
% Damping_Data.mass_factor = 10;
% Damping_Data.stiffness_factor = 0;
% 
% 
% Force_Data.type = "modal";
% Force_Data.mode_number = 1;
% Force_Data.continuation_variable = "frequency";
% Force_Data.amplitude = 2;
% % 
% 
% %--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 3e1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 10;
% Continuation_Opts.backward_steps = 500;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.collation_degree = 6;
% Continuation_Opts.parameter_range = [0.009,0.0209];
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% 
% 
% Dyn_Data = Dyn_Data.get_fe_output("forced_response",2,75);

% set_visualisation_level(3)
% Static_Data_S = load_static_data("clamped_beam_1");
% Static_Data_P = load_static_data("clamped_beam_copy_1");
% Static_Data_PF = load_static_data("clamped_beam_fixed_p_1");
% 
% Static_Data_S = Static_Data_S.add_validation_data(1:10);
% Static_Data_P = Static_Data_P.add_validation_data(1:10);
% Static_Data_PF = Static_Data_PF.add_validation_data(1:10);
% 
% plot_span = 1:4;
% span_length = size(plot_span,2);
% plot_outputs = zeros(span_length^2,2);
% plot_counter = 0;
% for iRow = 1:span_length
%     for iCol = 1:span_length
%         plot_counter = plot_counter + 1;
%         plot_outputs(plot_counter,:) = plot_span([iRow,iCol]);
%     end
% end
% 
% ax = plot_static_data("h_stiffness",Static_Data_PF,"outputs",plot_outputs);
% plot_static_data("h_stiffness",Static_Data_S,"outputs",plot_outputs,"axes",ax);