clear
% close all
set_visualisation_level(1)
set_logging_level(2)

system_name = "mems_arch_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = "all";
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% 
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e0;
Continuation_Opts.max_inc = 1e0;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);


% compare_validation(Dyn_Data,"validation error",2,"all");
% Dyn_Data = Dyn_Data.validate_solution(1,[6,11]);


Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 5e-2;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 200;
Dyn_Data = Dyn_Data.add_orbits(1,[20,22],"opts",Continuation_Opts);


% compare_solutions("energy",Dyn_Data,1);
% compare_validation(Dyn_Data,"validation error",1,"all");
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
quality_factor = 500;
omega_I = sqrt(Dyn_Data.Dynamic_Model.Model.reduced_eigenvalues(1));

Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = omega_I/quality_factor;
Damping_Data.stiffness_factor = 0;


Force_Data.type = "modal";
Force_Data.mode_number = 1;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


Force_Data.continuation_variable = "amplitude";
Force_Data.force_points = [0.0025,0.005,0.0075,0.01]*1e6;
% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e4;
Continuation_Opts.max_inc = 1e5;
Continuation_Opts.min_inc = 1e3;
Continuation_Opts.forward_steps = 100;
Continuation_Opts.backward_steps = 0;
%-----------------------------------------%
frequency_range = [2.6586e+06,2.7540e+06];
Force_Data.frequency = frequency_range(1);
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
Force_Data.frequency = frequency_range(2);
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
% %-------------------------------------------------------------------------%
% %-------------------------------------------------------------------------%

period_range = 2*pi*[0.5,1.5]./frequency_range;

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 1e3;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 174;
Continuation_Opts.parameter_range = period_range;
%-----------------------------------------%
Dyn_Data = Dyn_Data.restart_point(2,2,"IC","opts",Continuation_Opts);


%----------------------------------------------------------------------------%

% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 1;
Continuation_Opts.min_inc = 1e-1;
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 100;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%
Dyn_Data = Dyn_Data.restart_point(4,96,"frf_to_bb","opts",Continuation_Opts);

%----------------------------------------------------------------------------%
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 100;
Continuation_Opts.backward_steps = 00;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
Continuation_Opts.energy_limit_multiplier = 1;
Continuation_Opts.parameter_range = frequency_range;
%-----------------------------------------%
Dyn_Data = Dyn_Data.restart_point(5,1,"po","opts",Continuation_Opts);