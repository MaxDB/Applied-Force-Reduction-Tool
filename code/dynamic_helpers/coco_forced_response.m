function coco_forced_response(t0,z0,p0,Rom,type,Continuation_Settings,solution_number,Additional_Output,Nonconservative_Input)
%Calculates frfs using coco.
ODE_TOLERACE = 1e-9;
RESONANT_PHASE = pi/2;

% Set increment options
h0 = Continuation_Settings.initial_inc;
h_min = Continuation_Settings.min_inc;
h_max = Continuation_Settings.max_inc;
inc_forward = Continuation_Settings.forward_steps;
inc_backward = Continuation_Settings.backward_steps;

%Define EOM
switch type
    case "rom"
        Eom_Input = Rom.get_solver_inputs("coco_frf",Nonconservative_Input);

        switch Nonconservative_Input.continuation_variable
            case "amplitude"
                T = Eom_Input.Applied_Force_Data.period;
                funcs = {@(t,z,amp) coco_forced_eom(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data),...
                    @(t,z,amp) coco_forced_eom_dx(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data),...
                    @(t,z,amp) coco_forced_eom_dA(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data),...
                    @(t,z,amp) coco_forced_eom_dt(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data)};


                coll_args = [funcs, {t0',z0', {'F'}, p0}];
                cont_args = { 1, {'F'}, [0,max(Nonconservative_Input.force_points)*1.05]};
            case "frequency"
                amp = Eom_Input.Applied_Force_Data.amplitude;
                %probably parameter jacobian needs to be tensor?
                funcs = {@(t,z,T) coco_forced_eom(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data),...
                    @(t,z,T) coco_forced_eom_dx(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data),...
                    @(t,z,T) coco_forced_eom_dTper(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data),...
                    @(t,z,T) coco_forced_eom_dt(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data)};
            
                T0 = 2*pi/p0;
                coll_args = [funcs, {t0',z0', {'T'}, T0}];
                cont_args = { 1, {'po.period', 'T'}, Continuation_Settings.parameter_range'};
        end
    case "fom"
        % Model = Rom.Model;
        % load("geometry\" + Model.system_name+ "\" + Model.system_name + ".mat","analytic_eom");
        % Eom_Input = analytic_eom.get_solver_inputs("free");
        %
        % funcs = {@(z,zeta) direct_eom(0,z,zeta,Eom_Input.modal_restoring_force),...
        %     @(z,zeta) direct_eom_dx(0,z,zeta,Eom_Input.modal_stiffness),...
        %     @(z,zeta) direct_eom_dzeta(0,z)};
end


% Continuation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
%Collation Settings

prob = coco_set(prob, 'coll', 'NTST',   Continuation_Settings.initial_discretisation_num);  % [10] %initial number of discretisation intervals
prob = coco_set(prob, 'coll', 'NCOL',   Continuation_Settings.collation_degree);            % [4] %degree of interpolating polynomial
prob = coco_set(prob, 'coll', 'NTSTMN',   Continuation_Settings.min_discretisation_num);    % [5] %min number of discretisation intervals
prob = coco_set(prob, 'coll', 'NTSTMX',   Continuation_Settings.max_discretisation_num);    % [100] %max number of discretisation intervals



%ODE Settings
prob = coco_set(prob, 'ode', 'RelTol', ODE_TOLERACE);
prob = coco_set(prob, 'ode', 'AbsTol', ODE_TOLERACE*1e-2);


prob = ode_isol2po(prob, '', coll_args{:});

if Nonconservative_Input.continuation_variable == "frequency"
    % Constrain period
    [data,uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
end

% Monitor Energy
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data = coco_energy_init_data(fdata, 'po.orb');
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);

switch type
    case "rom"
        energy_func = @(prob,data,u) coco_energy(prob,data,u,Eom_Input.Potential_Polynomial,Eom_Input.Disp_Data,Eom_Input.input_order,"forced");
    case "fom"
        energy_func = @(prob,data,u) direct_energy(prob,data,u,Eom_Input);
end
prob = coco_add_func(prob, 'energy_monitor', energy_func, data, 'regular', 'ENERGY', 'uidx', uidx,'remesh',@coco_energy_remesh);

energy_limit = Rom.Model.energy_limit*Continuation_Settings.energy_limit_multiplier;
prob = coco_add_event(prob, 'VBP', 'boundary','ENERGY',energy_limit);
prob = coco_add_slot(prob, 'energy_slot',@coco_energy_print,data,'cont_print');

% Monitor additional output
switch Additional_Output.type
    case "physical displacement"
        node_map = Rom.Model.node_mapping;

        dof = node_map(Additional_Output.dof,2);
        Disp_Poly = Rom.Physical_Displacement_Polynomial;

        displacement_func = @(prob,data,u) coco_displacement(prob,data,u,Disp_Poly,dof);

        prob = coco_add_func(prob, 'displacement_monitor', displacement_func, data, 'regular', 'DISP', 'uidx', uidx,'remesh',@coco_energy_remesh);
        disp_points = Additional_Output.special_points;
        prob = coco_add_event(prob, 'X', 'special point','DISP',disp_points);
end

switch Nonconservative_Input.continuation_variable
    case "amplitude"
        prob = coco_add_event(prob, 'IC', 'F', Nonconservative_Input.force_points);
    case "frequency"
        

        % Detect resonance via phase difference
        phase_func = @(prob,data,u) coco_phase(prob,data,u,Eom_Input.Applied_Force_Data);
        prob = coco_add_func(prob, 'phase_monitor', phase_func, data, 'regular', 'PHASE', 'uidx', uidx,'remesh',@coco_energy_remesh);

        prob = coco_add_event(prob, 'RES', 'special point','PHASE',RESONANT_PHASE);
end


%Corrector Settings
prob = coco_set(prob,'corr','SubItMX', 10); % [4] number of damping steps
prob = coco_set(prob,'corr','ItMX', 50);    % [10] Maximum number of retries for an iteration 

%Continuation settings
prob = coco_set(prob, 'cont', 'NAdapt', 1);
prob = coco_set(prob, 'cont', 'NPR',    1);         % Number of steps between prints
prob = coco_set(prob, 'cont', 'ItMX',   [inc_backward,inc_forward]);  	% Number of continuation steps [backwards,forwards]

prob = coco_set(prob, 'cont', 'h0',     h0);
prob = coco_set(prob, 'cont', 'h_min',  h_min);
prob = coco_set(prob, 'cont', 'h_max',  h_max);

file_name = "temp\dynamic_sol_" + solution_number;

coco(prob, file_name, [], cont_args{:});
end