clear
close all
set_visualisation_level(1)
set_logging_level(2)

system_name = "exhaust_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
% Additional_Output.output = "physical displacement";
% Additional_Output.type = "max";
% Additional_Output.dof = 1563;
% Additional_Output.special_points = [0.25,0.5,0.75,1,1.5,2,3]*1.5e-3;
% Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e1;
% Continuation_Opts.max_inc = 1e1;
% Continuation_Opts.min_inc = 1e1;
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 6;
%-----------------------------------------%

% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% compare_validation(Dyn_Data,"validation error",1,"all")
% Dyn_Data = Dyn_Data.validate_solution(1,"all");
% Dyn_Data = Dyn_Data.get_fe_output("periodicity",1,"X");
% Dyn_Data = Dyn_Data.get_max_disp_stress(1,46);


%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.frequency_range = [800,1300];
%-----------------------------------------%


Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 0;
Damping_Data.stiffness_factor = 4e-5;

Force_Data.type = "point force";
Force_Data.dof = 1563;
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 850;
Force_Data.amplitude = 100;

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom");
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom","method","fc");
