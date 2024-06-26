function [r,theta,f,E,additional_data] = add_point_matlab(applied_force,add_data_type,Model)
num_r_modes = length(Model.reduced_modes);
%load system
Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
static_equation = Analytic_Eom.get_static_equation;
potential_equation = Analytic_Eom.get_potential_equation;
linear_stiffness = Model.stiffness;

%define force
force_transform = Model.mass*Model.reduced_eigenvectors;

%preallocation
num_loadcases = size(applied_force,2);

f = zeros(num_r_modes,num_loadcases);
displacement = zeros(Model.num_dof,num_loadcases);
E = zeros(1,num_loadcases);
switch add_data_type
    case "none"
        additional_data = [];
    case "stiffness"
        additional_data = zeros(Model.num_dof,Model.num_dof,num_loadcases);
    case "perturbation"
        error("Not implemented yet for direct systems")
    otherwise
        error("additional data type not recognised")
end

%solver settings
options = optimoptions('fsolve','Display','off');


for iLoad = 1:num_loadcases
    step_force = applied_force(:,iLoad);

    physical_force = force_transform*step_force;
    disp_guess = linear_stiffness\physical_force;

    [physical_disp,~,~,~,stiffness] = fsolve(@(x) static_equation(x) - physical_force,disp_guess,options);

    potential = potential_equation(physical_disp);

    f(:,iLoad) = step_force;
    displacement(:,iLoad) = physical_disp;
    E(1,iLoad) = potential;

    switch add_data_type
        case "none"

        case "stiffness"
            additional_data(:,:,iLoad) = stiffness; %#ok<AGROW>
    end

end

%analyse data

disp_transform = force_transform';
r = disp_transform*displacement;
% theta = displacement - Model.reduced_eigenvectors*r;
theta = displacement;
end
