clear
close all
set_visualisation_level(1)

system_name = "clamped_beam_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 248;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

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

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

%-------------------
% %Low damping
% target_damping = 0.01;
% natural_freq = sqrt(Dyn_Data.Dynamic_Model.Model.reduced_eigenvalues);
% 
% Damping_Data.damping_type = "rayleigh";
% Damping_Data.mass_factor = 0;
% Damping_Data.stiffness_factor = 2*target_damping/natural_freq(1);
% 
% Force_Data.type = "point";
% Force_Data.dof = 248;
% Force_Data.continuation_variable = "frequency";
% Force_Data.frequency = 350;
% Force_Data.amplitude = 0.05;
% 
% % --------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 2e-2;
% Continuation_Opts.max_inc = 2e-2;
% Continuation_Opts.min_inc = 5e-3;
% 
% Continuation_Opts.forward_steps = 200;
% Continuation_Opts.backward_steps = 200;
% Continuation_Opts.collocation_degree = 6;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.frequency_range = [290,450];
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom");

%------
% Force_Data.type = "modal";
% Force_Data.continuation_variable = "frequency";
% Force_Data.mode_number = 1;
% Force_Data.frequency = 350;
% Force_Data.amplitude = 0.2488;
% 
% % --------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 2e-2;
% Continuation_Opts.max_inc = 2e-2;
% Continuation_Opts.min_inc = 5e-3;
% 
% Continuation_Opts.forward_steps = 200;
% Continuation_Opts.backward_steps = 200;
% Continuation_Opts.collocation_degree = 6;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.frequency_range = [290,450];
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom","method","fc");


%-------------------
%point force

% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 2e-2;
Continuation_Opts.max_inc = 2e-2;
Continuation_Opts.min_inc = 5e-3;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.frequency_range = [290,450];
%-----------------------------------------%

Force_Data.type = "point";
Force_Data.dof = 248;
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 350;
Damping_Data.damping_type = "rayleigh";

%2,43
%--
Force_Data.amplitude = 0.05;
target_damping = 0.01;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,0]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = 0;

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%--

% Force_Data.amplitude = 1;
% target_damping = 0.2;
% damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
% Damping_Data.mass_factor = damping_coeffs(1);
% Damping_Data.stiffness_factor = damping_coeffs(2);
% 
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%--
% Force_Data.amplitude = 5;
% target_damping = 1;
% damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
% Damping_Data.mass_factor = damping_coeffs(1);
% Damping_Data.stiffness_factor = damping_coeffs(2);
% 
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);

%------
% Force_Data.type = "modal";
% Force_Data.continuation_variable = "frequency";
% Force_Data.mode_number = 1;
% Force_Data.frequency = 350;
% Force_Data.amplitude = 1;
% 
% % --------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 2e-2;
% Continuation_Opts.max_inc = 2e-2;
% Continuation_Opts.min_inc = 5e-3;
% 
% Continuation_Opts.forward_steps = 200;
% Continuation_Opts.backward_steps = 200;
% Continuation_Opts.collocation_degree = 6;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.frequency_range = [290,450];
% %-----------------------------------------%
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom","method","fc");
% 
% 


%--------
% %Validation
% Dyn_Data = initalise_dynamic_data("clamped_beam_11001");
% Dyn_Data = Dyn_Data.validate_solution(2,3);