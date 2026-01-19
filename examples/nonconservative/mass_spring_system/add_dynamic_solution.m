clear
% close all
set_visualisation_level(1)

system_name = "mass_spring_roller_1";
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 2;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;

Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.frequency_points = [1.86,3];


%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts,"type","fom");

%-----------------------------------
Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 0.1;
Damping_Data.stiffness_factor = 1e-7;


Force_Data.type = "point force";
Force_Data.dof = 3;
Force_Data.continuation_variable = "frequency";
Force_Data.amplitude = 0.5;
Force_Data.frequency = 1.1;


% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 1e-3;

Continuation_Opts.backward_steps = 500;
Continuation_Opts.frequency_range = [0.995,3.2];
%-----------------------------------------%
% 
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","fom");