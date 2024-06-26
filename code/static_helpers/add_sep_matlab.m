function [r,theta,f,E,additional_data,sep_id] = ...
    add_sep_matlab(force_ratio,num_loadcases,add_data_type,Model)

num_seps = size(force_ratio,2);
num_r_modes = length(Model.reduced_modes);
%load system
Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
static_equation = Analytic_Eom.get_static_equation;
potential_equation = Analytic_Eom.get_potential_equation;
linear_stiffness = Model.stiffness;

%define force
applied_force = force_ratio/num_loadcases;
force_transform = Model.mass*Model.reduced_eigenvectors;

%preallocation
total_static_steps = 2*num_seps*num_loadcases;

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

%solver settings
options = optimoptions('fsolve','Display','off');

load_step_counter = 0;

for iSep = 1:num_seps
        sep_force = applied_force(:,iSep);

        physical_force = force_transform*sep_force;
        disp_guess = linear_stiffness\physical_force;

        reached_energy_limit = 0;
        step_force = zeros(size(physical_force));
        modal_step_force = zeros(num_r_modes,1);
        while ~reached_energy_limit
            load_step_counter = load_step_counter + 1;
            step_force = step_force + physical_force;
            modal_step_force = modal_step_force + sep_force;

            [physical_disp,~,~,~,stiffness] = fsolve(@(x) static_equation(x) - step_force,disp_guess,options);
            disp_guess = physical_disp*(load_step_counter+1)/load_step_counter;

            potential = potential_equation(physical_disp);
            
            reached_energy_limit = potential > Model.fitting_energy_limit;

            sep_id(1,load_step_counter) = iSep;
            f(:,load_step_counter) = modal_step_force;
            displacement(:,load_step_counter) = physical_disp;
            E(1,load_step_counter) = potential;

            switch add_data_type
                case "none"

                case "stiffness"
                    additional_data(:,:,load_step_counter) = stiffness; %#ok<AGROW>
            end
        end
end

%analyse data
remove_index = (load_step_counter+1):size(sep_id,2);
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