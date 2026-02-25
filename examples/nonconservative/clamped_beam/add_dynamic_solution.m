clear
close all
set_visualisation_level(1)

system_name = "clamped_beam_11001";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Additional_Output.output = "physical displacement";
% Additional_Output.type = "max";
% Additional_Output.dof = 182;
% Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
%%% one mode
% Continuation_Opts.initial_inc = 5e-2;
% Continuation_Opts.max_inc = 5e-2;
% Continuation_Opts.min_inc = 5e-3;
%%% two mode
Continuation_Opts.initial_inc = 2e-2;
Continuation_Opts.max_inc = 2e-2;
Continuation_Opts.min_inc = 5e-3;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.min_discretisation_num = 20;
%-----------------------------------------%
% 
% Dyn_Data = Dyn_Data.add_backbone(1001,"opts",Continuation_Opts);
% % compare_validation(Dyn_Data,"validation error",1,"all")
% Dyn_Data.validate_solution(1,2)
% 
% %-------------------------------------------
% Continuation_Opts.initial_inc = 2e-3;
% Continuation_Opts.min_inc = 1e-3;
% Continuation_Opts.max_inc = 2e-3;
% Continuation_Opts.forward_steps = 0;
% Continuation_Opts.backward_steps = 250;
% %-----------------------------------------
% Dyn_Data = Dyn_Data.add_orbits(1,[7,11],"opts",Continuation_Opts);


Damping_Data.damping_type = "rayleigh";
% Damping_Data.mass_factor = 16;
Damping_Data.mass_factor = 0;
Damping_Data.stiffness_factor = 1e-3;


% Force_Data.type = "point force";
% Force_Data.dof = 362;
% Force_Data.continuation_variable = "frequency";
% Force_Data.frequency = 350;
% Force_Data.amplitude = 0.4;

% Force_Data.type = "point force";
% Force_Data.dof = 182;
% Force_Data.continuation_variable = "frequency";
% Force_Data.frequency = 350;
% Force_Data.amplitude = 0.5;


Force_Data.type = "uniform";
Force_Data.direction = 2;
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 350;
Force_Data.amplitude = 0.5;

% Force_Data.type = "modal";
% Force_Data.mode_number = 1;
% Force_Data.continuation_variable = "frequency";
% Force_Data.frequency = 350;
% Force_Data.amplitude = 3;


% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 2e-2;
Continuation_Opts.max_inc = 2e-2;
Continuation_Opts.min_inc = 5e-3;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.frequency_range = [290,600];
%-----------------------------------------%
% 
%% RE ENABLE JACOBIAN ETC.
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom");
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom","method","fc");


Dyn_Data = Dyn_Data.get_fe_output("forced_response",1,54);
