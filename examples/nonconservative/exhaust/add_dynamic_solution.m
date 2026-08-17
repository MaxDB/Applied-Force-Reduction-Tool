clear
close all
set_visualisation_level(1)
set_logging_level(2)

system_name = "exhaust_1";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 1563;
Additional_Output.special_points = [0.25,0.5,0.75,1,1.5,2,3]*1.5e-3;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;

Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.min_discretisation_num = 20;
%-----------------------------------------%

if size(Dyn_Data,1) == 0
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
end
% compare_validation("exhaust_1","energy",1,"all");

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 500;
Continuation_Opts.frequency_range = [1100,1400];
Continuation_Opts.collocation_degree = 8;
%-----------------------------------------%

Force_Data.type = "uniform";
Force_Data.direction = 3;
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 1110;

Damping_Data.damping_type = "rayleigh";

modal_damping = 0.000425;
Damping_Data.mass_factor = 2*sqrt(Dyn_Data.Dynamic_Model.Model.reduced_eigenvalues)*modal_damping; 
Damping_Data.stiffness_factor = 0;

%-
Force_Data.amplitude = 0.5;

Force_Data.frequency = 1110;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);

Force_Data.frequency = 1300;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);
%-
Force_Data.amplitude = 1;


Force_Data.frequency = 1110;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);

Force_Data.frequency = 1300;
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);


