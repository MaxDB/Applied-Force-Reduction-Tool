clear
close all
set_visualisation_level(1)

system_name = "clamped_beam_131001";
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
if size(Dyn_Data,1) == 0 && all(Dyn_Data.Dynamic_Model.Model.reduced_modes < 1000)
Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
end
%-------------------
%-------------------
%point force

% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 2e-2;
Continuation_Opts.max_inc = 2e-2;
Continuation_Opts.min_inc = 1e-3;

Continuation_Opts.forward_steps = 1000;
Continuation_Opts.backward_steps = 1000;
Continuation_Opts.collocation_degree = 8;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.frequency_range = [290,600];
%-----------------------------------------%

%-----------------------------------------%
Force_Data.type = "point";
Force_Data.dof = 248;
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 500;
Damping_Data.damping_type = "rayleigh";
%-----------------------------------------%

%--
Force_Data.amplitude = 5;
target_damping = 0.3;
%-
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,0]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = 0;

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);


