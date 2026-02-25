clear
% close all
set_visualisation_level(1)

system_name = "mass_spring_axial_1";
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 10;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e1;
Continuation_Opts.max_inc = 1e2;
Continuation_Opts.min_inc = 1e0;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;

Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 8;

%-----------------------------------------%
% 
% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts,"type","fom");
% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts,"type","rom");

%-----------------------------------
% Damping_Data.damping_type = "rayleigh";
% Damping_Data.mass_factor = 0;
% Damping_Data.stiffness_factor = 0.0035;

damping_coeff = 100;
Damping_Data.damping_type = "physical";
Damping_Data.damping_matrix = get_chain_damping_matrix(damping_coeff,Dyn_Data.Dynamic_Model.Model);


% Force_Data.type = "shape";
% Force_Data.shape = ones(10,1);
% Force_Data.continuation_variable = "frequency";
% Force_Data.frequency = 20;
% Force_Data.amplitude = 200;

Force_Data.type = "point";
Force_Data.dof = 10;
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 24;
Force_Data.amplitude = 1500;




% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e1;
Continuation_Opts.max_inc = 1e2;
Continuation_Opts.min_inc = 1e0;

Continuation_Opts.backward_steps = 500;
Continuation_Opts.frequency_range = [5,25];
%-----------------------------------------%
% 

% % 
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom");
Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","fom");




function damping_matrix = get_chain_damping_matrix(params,Model)
num_dof = Model.num_dof;
diag_matrix = zeros(num_dof+2);
diag_index = ((1:num_dof)'.*[1,1])+1;
diag_right_index = diag_index + [0,1];
diag_left_index = diag_index + [0,-1];

for iDof = 1:num_dof
    diag_matrix(diag_index(iDof,1),diag_index(iDof,2)) = 2;
    diag_matrix(diag_left_index(iDof,1),diag_left_index(iDof,2)) = -1;
    diag_matrix(diag_right_index(iDof,1),diag_right_index(iDof,2)) = -1;
end
diag_matrix(1,:) = [];
diag_matrix(:,1) = [];
diag_matrix(end,:) = [];
diag_matrix(:,end) = [];

diag_matrix(end,end) = 1; % free BC

damping_matrix = params*diag_matrix;
end