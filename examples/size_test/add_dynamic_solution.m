clear

close all
set_visualisation_level(1)
set_logging_level(3)

system_name = "mems_arch_16";
Dyn_Data = initalise_dynamic_data(system_name);
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [0,36e-3,10e-3];
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 5e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 2500;
Continuation_Opts.backward_steps = 0;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
% -----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);



%---
potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2.69e6,[1.5e-7,1e-7]);

Continuation_Opts.collation_degree = 10;
Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);

%--
% Continuation_Opts.initial_inc = 5e-2;
% Continuation_Opts.max_inc = 5e-2;
% Continuation_Opts.min_inc = 5e-2;
% Continuation_Opts.forward_steps = 500;
% Continuation_Opts.backward_steps = 00;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.collation_degree = 6;
% %--
% Dyn_Data = Dyn_Data.add_orbits(2,[5,8],"opts",Continuation_Opts);
% 
% 
% 
% %----
% %high denisty backbone
% % --------- Continuation Settings ---------%
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 1e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 2500;
% Continuation_Opts.backward_steps = 0;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.collation_degree = 6;
% % -----------------------------------------%
% 
% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
% 
% Continuation_Opts.backward_steps = 2;
% Dyn_Data = Dyn_Data.restart_point(2,3,"po","opts",Continuation_Opts);
% 
% 
% % -----------------------------------------%
% Dyn_Data_16 = initalise_dynamic_data("mems_arch_16");
% Orbit = Dyn_Data_16.get_orbit(2,1);
% t = Orbit.tbp';
% r_16 = Orbit.xbp(:,1:2)';
% r_dot_16 = Orbit.xbp(:,3:4)';
% 
% Rom_16 = Dyn_Data_16.Dynamic_Model;
% x = Rom_16.expand(r_16);
% x_dot = Rom_16.expand_velocity(r_16,r_dot_16);
% 
% r_evec = Dyn_Data.Dynamic_Model.Model.reduced_eigenvectors;
% mass = Dyn_Data.Dynamic_Model.Model.mass;
% transform = r_evec'*mass;
% 
% r = transform*x;
% r_dot = transform*x_dot;
% 
% orbit_ic = {t,[r;r_dot]};
% 
% % -----------------------------------------%
% Continuation_Opts.initial_inc = 5e-1;
% Continuation_Opts.max_inc = 5e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 0;
% 
% Dyn_Data = Dyn_Data.add_backbone(1,"ic",orbit_ic,"opts",Continuation_Opts);
% -----------------------------------------%

% Continuation_Opts.initial_inc = 1e0;
% Continuation_Opts.max_inc = 1e0;
% Continuation_Opts.min_inc = 1e0;

% Continuation_Opts.forward_steps = 10;
% Continuation_Opts.backward_steps = 0;

% Continuation_Opts.initial_inc = 2e0;
% Continuation_Opts.max_inc = 2e0;
% Continuation_Opts.min_inc = 2e0;
% 
% Continuation_Opts.forward_steps = 0;
% Continuation_Opts.backward_steps = 10;
% Dyn_Data = Dyn_Data.restart_point(2,7,"po","opts",Continuation_Opts);
% % -----------------------------------------%
% Continuation_Opts.initial_inc = 5e-1;
% Continuation_Opts.max_inc = 5e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 100;
% Continuation_Opts.backward_steps = 100;
% % 
% % Dyn_Data = Dyn_Data.restart_point(3,4,"po","opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.restart_point(3,2,"po","opts",Continuation_Opts);
% Dyn_Data = Dyn_Data.remove_solution(3);
% %-----
% 
% Continuation_Opts.initial_inc = 1e-1;
% Continuation_Opts.max_inc = 1e-1;
% Continuation_Opts.min_inc = 1e-2;
% Continuation_Opts.forward_steps = 0;
% Continuation_Opts.backward_steps = 50;
% Continuation_Opts.initial_discretisation_num = 20;
% Continuation_Opts.max_discretisation_num = 250;
% Continuation_Opts.min_discretisation_num = 20;
% Continuation_Opts.collation_degree = 6;
% %---
% Dyn_Data = Dyn_Data.add_orbits(1,[4,6],"opts",Continuation_Opts);
% 
% Continuation_Opts.forward_steps = 50;
% Continuation_Opts.backward_steps = 0;
% Dyn_Data = Dyn_Data.add_orbits(2,[11,13],"opts",Continuation_Opts);