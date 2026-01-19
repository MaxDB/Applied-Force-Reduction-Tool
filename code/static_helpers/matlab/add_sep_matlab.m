function [reduced_displacement,physical_displacement,f,E,additional_data,sep_id] = ...
    add_sep_matlab(force_ratio,num_loadcases,add_data_type,Model)


conservative = true;
if iscell(Model)
    [Model,Nc_Data] = Model{:};
   

    conservative = isempty(Nc_Data);
end

[num_modes,num_seps] = size(force_ratio);

%load system
Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
static_equation = Analytic_Eom.get_static_equation;
potential_equation = Analytic_Eom.get_potential_equation;
linear_stiffness = Model.stiffness;

%define force
applied_force = force_ratio/num_loadcases;
force_transform = Model.mass*Model.reduced_eigenvectors;

if ~conservative
    nc_force_transform = Nc_Data.max_amplitude.*Nc_Data.force_shape;
    force_transform = [force_transform,nc_force_transform];
    num_r_modes = size(Model.reduced_modes,1);
    num_applied_force = Nc_Data.num_applied_forces;
end

%preallocation
total_static_steps = 2*num_seps*num_loadcases;

sep_id = zeros (1,total_static_steps);
f = zeros(num_modes,total_static_steps);
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
fsolve_options = optimoptions('fsolve','Display','off');

load_step_counter = 0;
nc_removal_indicies =[];
for iSep = 1:num_seps
        sep_force = applied_force(:,iSep);

        physical_force = force_transform*sep_force;
        disp_guess = linear_stiffness\physical_force;
        disp_guess(isnan(disp_guess)) = 0;

        reached_limit = 0;
        step_force = zeros(size(physical_force));
        modal_step_force = zeros(num_modes,1);
        while ~reached_limit
            load_step_counter = load_step_counter + 1;
            step_force = step_force + physical_force;
            modal_step_force = modal_step_force + sep_force;

            [physical_disp,~,~,~,stiffness] = fsolve(@(x) static_equation(x) - step_force,disp_guess,fsolve_options);
            disp_guess = physical_disp*(load_step_counter+1)/load_step_counter;

            potential = potential_equation(physical_disp);
            
            reached_limit = potential > Model.fitting_energy_limit;

            if ~conservative
                force_amplitudes = modal_step_force((num_r_modes+1):end);
                if  any(abs(force_amplitudes)>1)
                    reached_limit = 1;
                    nc_removal_indicies = [nc_removal_indicies,load_step_counter]; %#ok<AGROW>
                end
                
            end

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
remove_index = [nc_removal_indicies,remove_index];
sep_id(:,remove_index) = [];
f(:,remove_index) = [];
displacement(:,remove_index) = [];
E(:,remove_index) = [];
switch add_data_type
    case "none"

    case "stiffness"
        additional_data(:,:,remove_index) = [];
end


r_transform = force_transform';
r = r_transform*displacement;
physical_displacement = displacement;
reduced_displacement = r;

if conservative
    return
end




sin_fraction = f((num_r_modes+1):end,:);

p = asin(sin_fraction);
x_p = nc_force_transform*sin(p);
x_r = physical_displacement - x_p;
r = r_transform*x_r;

r((num_r_modes+1):end,:) = p;
reduced_displacement = r;

% reduced_displacement = [r;p];

% reduced_displacement = r;

% r = r(1:num_r_modes,:);
% x_r = Model.reduced_eigenvectors*r;
% x_p = displacement-x_r;
% p = nc_force_transform'*x_p;
% 
% reduced_displacement = [r;p];
end