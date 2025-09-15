function [r,theta,f,E,additional_data,sep_id] = ...
    add_sep_matlab_riks(force_ratio,target_loadcases,add_data_type,Model)

num_seps = size(force_ratio,2);
num_r_modes = length(Model.reduced_modes);
%load system
Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
static_equation = Analytic_Eom.get_static_equation;
potential_equation = Analytic_Eom.get_potential_equation;
stiffness_equation = Analytic_Eom.get_stiffness;


System_Eqs.static = static_equation;
System_Eqs.potential = potential_equation;
System_Eqs.stiffness = stiffness_equation;

energy_limit = Model.fitting_energy_limit;
%define force
applied_force = force_ratio/target_loadcases;
force_transform = Model.mass*Model.reduced_eigenvectors;

%preallocation
total_static_steps = 2*num_seps*target_loadcases;

sep_id = zeros (1,total_static_steps);
f = zeros(num_r_modes,total_static_steps);
displacement = zeros(Model.num_dof,total_static_steps);
E = zeros(1,total_static_steps);
switch add_data_type
    case "none"
        additional_data = [];
    case "stiffness"
        additional_data = zeros(Model.num_dof,Model.num_dof,total_static_steps);
    case "perturbation"
        error("Not implemented yet for direct systems")
    otherwise
        error("additional data type not recognised")
end

% linear_stiffness = Model.stiffness;
% fsolve_options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
%solver settings
load_step_end = 0;
for iSep = 1:num_seps
        sep_force = applied_force(:,iSep);

        physical_force = force_transform*sep_force;
        % disp_guess = linear_stiffness\physical_force;
        % disp_guess(isnan(disp_guess)) = 0;
        % [lin_lambda] = fsolve(@(lambda) static_equation(lambda*disp_guess) - lambda*physical_force,1,fsolve_options);
        % initial_point = {lin_lambda*disp_guess,lin_lambda};
        target_force = physical_force*target_loadcases;

        initial_condition = {0,0};
        limit_reached = 0;
        while ~limit_reached

            [disp_sep,lambda_sep,potential_sep,stiffness_sep] = find_sep_system(System_Eqs,target_force,target_loadcases,energy_limit,"ic",initial_condition);
            num_sep_steps = size(lambda_sep,2);
            load_step_start = load_step_end + 1;
            load_step_end = load_step_end + num_sep_steps;
            load_span = load_step_start:load_step_end;

            modal_step_force = lambda_sep.*sep_force*target_loadcases;

            sep_id(1,load_span) = iSep;
            f(:,load_span) = modal_step_force;
            displacement(:,load_span) = disp_sep;
            E(1,load_span) = potential_sep;


            switch add_data_type
                case "none"

                case "stiffness"
                    additional_data(:,:,load_span) = stiffness_sep; %#ok<AGROW>
            end

            limit_reached = potential_sep(end) > Model.energy_limit;
            initial_condition = {disp_sep(:,end),lambda_sep(end)};
        end

end

%analyse data
excess_index = (load_step_end+1):size(sep_id,2);
limit_index = E > energy_limit;
remove_index = limit_index;
remove_index(excess_index) = true;

sep_id(:,remove_index) = [];
f(:,remove_index) = [];
displacement(:,remove_index) = [];
E(:,remove_index) = [];
switch add_data_type
    case "none"

    case "stiffness"
        additional_data(:,:,remove_index) = [];
end

disp_transform = force_transform';
r = disp_transform*displacement;
% theta = displacement - Model.reduced_eigenvectors*r;
theta = displacement;
end