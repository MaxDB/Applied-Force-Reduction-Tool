clear
close all
set_visualisation_level(1)
set_logging_level(3)

system_name = "cantilever_f_2001";
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 722;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 6;
%-----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
%-----------------------------------------%
% 
% Dyn_Data_25 = initalise_dynamic_data("cantilever_comp_1");
% Orbit = Dyn_Data_25.get_orbit(3,31);
% t0 = Orbit.tbp;
% z0 = Orbit.xbp';
% 
% Dyn_Data = Dyn_Data.add_backbone(1,"ic",{t0,z0},"opts",Continuation_Opts);


% Dyn_Data = Dyn_Data.validate_solution(1,2001);
% compare_validation(Dyn_Data,"validation error",1,{"1f"})


% Dyn_Data = Dyn_Data.get_fe_output("periodicity",1,42);
