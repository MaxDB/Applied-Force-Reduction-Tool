function coco_backbone(t0,z0,Rom,type,Continuation_Settings,solution_number,Additional_Output,varargin)
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

Restart_Data.initial_solution_type = "initial_solution";

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "restart"
            Restart_Data = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
%Calculates bb curves using coco.
ODE_TOLERACE = 1e-9;

% Set increment options
h0 = Continuation_Settings.initial_inc;
h_min = Continuation_Settings.min_inc;
h_max = Continuation_Settings.max_inc;
inc_forward = Continuation_Settings.forward_steps;
inc_backward = Continuation_Settings.backward_steps;

%Define EOM
switch type
    case "rom"
        Eom_Input = Rom.get_solver_inputs("coco_backbone");

        

        if Continuation_Settings.inertial_compensation
            funcs = {@(z,zeta) coco_eom(0,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data),...
                @(z,zeta) coco_eom_dx(0,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data),...
                @(z,zeta) coco_eom_dzeta(0,z)};
        else
            funcs = {@(z,zeta) ice_eom(0,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data),...
                @(z,zeta) ice_eom_dx(0,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data),...
                @(z,zeta) ice_eom_dzeta(0,z)};
        end

    case "fom"
        Model = Rom.Model;
        Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
        Eom_Input = Analytic_Eom.get_solver_inputs("free");

        funcs = {@(z,zeta) direct_eom(0,z,zeta,Eom_Input.modal_restoring_force),...
            @(z,zeta) direct_eom_dx(0,z,zeta,Eom_Input.modal_stiffness),...
            @(z,zeta) direct_eom_dzeta(0,z)};
end


% Continuation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob = coco_prob();
% prob = coco_set(prob, 'ode', 'autonomous', false);



prob = coco_set(prob, 'coll', 'NTST',   Continuation_Settings.initial_discretisation_num);  % [10] %initial number of discretisation intervals
prob = coco_set(prob, 'coll', 'NCOL',   Continuation_Settings.collocation_degree);            % [4] %degree of interpolating polynomial
prob = coco_set(prob, 'coll', 'NTSTMN',   Continuation_Settings.min_discretisation_num);    % [5] %min number of discretisation intervals
prob = coco_set(prob, 'coll', 'NTSTMX',   Continuation_Settings.max_discretisation_num);    % [100] %max number of discretisation intervals

% cont_args = { 1, {'po.period', 'zet', 'ENERGY'}, [0,inf]};
cont_args = { 1, {'po.period', 'zet'}, []};

%ODE Settings
prob = coco_set(prob, 'ode', 'RelTol', ODE_TOLERACE);
prob = coco_set(prob, 'ode', 'AbsTol', ODE_TOLERACE*1e-2);

switch Restart_Data.initial_solution_type
    case "initial_solution"
        prob = coco_set(prob, 'cont', 'NAdapt', 1);
        %collocation Settings
        coll_args = [funcs, {t0',z0', {'zet'}, 0}];
        prob = ode_isol2po(prob, '', coll_args{:});
    case "branch_point"
        prob = coco_set(prob, 'cont', 'NAdapt', 5);
        restarted_file_name = Rom.data_path + "dynamic_sol_" + Restart_Data.sol_num;
        prob = ode_BP2po(prob, '', convertStringsToChars(restarted_file_name), Restart_Data.orbit_id);
    case "period_doubling"
        prob = coco_set(prob, 'cont', 'NAdapt', 5);
        restarted_file_name = Rom.data_path + "dynamic_sol_" + Restart_Data.sol_num;
        prob = ode_PD2po(prob, '', convertStringsToChars(restarted_file_name), Restart_Data.orbit_id);
end

% Constrain period -> only necessary if set to non-autonomous
% [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
% maps = data.coll_seg.maps;
% prob = coco_add_pars(prob, 'section', uidx(maps.x0_idx(1)), 'y0');

% Monitor Energy
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data = coco_energy_init_data(fdata, 'po.orb');
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);

switch type
    case "rom"
        energy_func = @(prob,data,u) coco_energy(prob,data,u,Eom_Input.Potential_Polynomial,Eom_Input.Disp_Data,Eom_Input.input_order,"free");
    case "fom"
        energy_func = @(prob,data,u) direct_energy(prob,data,u,Eom_Input,1);
end
prob = coco_add_func(prob, 'energy_monitor', energy_func, data, 'regular', 'ENERGY', 'uidx', uidx,'remesh',@coco_energy_remesh);

energy_limit = Rom.Model.energy_limit*Continuation_Settings.energy_limit_multiplier;
prob = coco_add_event(prob, 'VBP', 'boundary','ENERGY',energy_limit);
prob = coco_add_slot(prob, 'energy_slot',@coco_energy_print,data,'cont_print');

% Monitor additional output
switch Additional_Output.output
    case "physical displacement"
        disp_func = Additional_Output.output_func;
        coco_disp_func = @(prob,data,u) coco_displacement(prob,data,u,disp_func);

        prob = coco_add_func(prob, 'displacement_monitor', coco_disp_func, data, 'regular', 'DISP', 'uidx', uidx,'remesh',@coco_energy_remesh);
        disp_points = Additional_Output.special_points;
        prob = coco_add_event(prob, 'X', 'special point','DISP',disp_points);
end

% Monitor frequency
% coco_frequency_func = @(prob,data,u) coco_frequency(prob,data,u);
% prob = coco_add_func(prob, 'frequency_monitor', coco_frequency_func, data, 'regular', 'FREQ', 'uidx', uidx,'remesh',@coco_energy_remesh);
freq_func = @(prob,data,u) coco_frequency(prob,data,u);
prob = coco_add_func(prob, 'frequency_monitor', freq_func, data, 'regular', 'FREQ', 'uidx', uidx,'remesh',@coco_energy_remesh);
prob = coco_add_event(prob, 'EP','boundary','FREQ',Continuation_Settings.parameter_range);

if ~isempty(Continuation_Settings.frequency_points)
    prob = coco_add_event(prob, 'X','special point','FREQ',Continuation_Settings.frequency_points);
end


%Corrector Settings
prob = coco_set(prob,'corr','SubItMX', 10); % [4] number of damping steps
prob = coco_set(prob,'corr','ItMX', 100);    % [10] Maximum number of retries for an iteration 

%Continuation settings

prob = coco_set(prob, 'cont', 'NPR',    1);         % Number of steps between prints
prob = coco_set(prob, 'cont', 'ItMX',   [inc_backward,inc_forward]);  	% Number of continuation steps [backwards,forwards]

prob = coco_set(prob, 'cont', 'h0',     h0);
prob = coco_set(prob, 'cont', 'h_min',  h_min);
prob = coco_set(prob, 'cont', 'h_max',  h_max);

file_name = "temp\dynamic_sol_" + solution_number;

coco(prob, file_name, [], cont_args{:});
end