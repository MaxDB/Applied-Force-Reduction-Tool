function Static_Data = sep_points_validation_new(Static_Data)
MAXIMUM_DEGREE = [11,11,9,7];
SEP_DENSITY_MAX = [inf,251,11,5];

Validation_Opts = Static_Data.Validation_Options;
max_interpolation_error = Validation_Opts.maximum_interpolation_error;
initial_degree = Validation_Opts.minimum_degree;
max_iterations = Validation_Opts.maximum_iterations;

Model = Static_Data.Model;
linear_mass = Model.mass;
Static_Opts = Model.Static_Options;
num_sep_loadcases = Validation_Opts.num_added_points;
max_sep_loadcases = Static_Opts.maximum_loadcases;
max_iteration_loadcases = Validation_Opts.max_added_points;



num_modes = length(Model.reduced_modes);
max_degree = MAXIMUM_DEGREE(num_modes);
degree_range = [initial_degree,max_degree];

solve_opts = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','none');

fitting_energy_limit = Model.fitting_energy_limit;
energy_limit = Model.energy_limit;

found_force_ratios = Static_Data.unit_sep_ratios;
num_original_seps = size(found_force_ratios,2);


sep_density = 1;
unit_force_ratios = found_force_ratios;

scaffold_points = Static_Data.scaffold_points;
data_available = any(scaffold_points == 0);
if data_available
    added_data = ~scaffold_points;
    [Static_Data,Static_Data_Added] = Static_Data.remove_loadcases(added_data);
end

for iIteration = 1:(max_iterations+1)

    validation_iteration_start = tic;

    sep_density = min(sep_density + 2,SEP_DENSITY_MAX(num_modes));
    new_unit_force_ratios = add_sep_ratios(num_modes,sep_density,found_force_ratios);
    unit_force_ratios = [unit_force_ratios,new_unit_force_ratios]; %#ok<AGROW>
    scaled_force_ratios = scale_sep_ratios(unit_force_ratios,Static_Data.Model.calibrated_forces);
    found_force_ratios = unit_force_ratios;

    num_validated_seps = size(unit_force_ratios,2);
    
    validated_seps = (1:num_validated_seps);
    
    

    V = Static_Data.potential_energy;
    validated_points = V < energy_limit;

    r = Static_Data.reduced_displacement(:,validated_points);
    f = Static_Data.restoring_force(:,validated_points);
    displacement = Static_Data.condensed_displacement(:,validated_points);


    [force_degree,worst_force_index] = minimum_polynomial_degree(Static_Data,"Force",r,f,degree_range); %#ok<ASGLU>
    [disp_degree,worst_disp_index] = minimum_polynomial_degree(Static_Data,"Condensed_Displacement",r,displacement,degree_range);
    
    if num_sep_loadcases == 0
        continue
    end
    %find predicted seps and compare
    check_interpolation_start = tic;

    Rom_One = Reduced_System(Static_Data,[force_degree,disp_degree],"full");
    validation_plot(1,Static_Data,Rom_One,worst_disp_index)
    Rom_Two = Reduced_System(Static_Data,[force_degree,disp_degree] + 2,"full");
    
    Potential_Poly_One = Rom_One.Potential_Polynomial;
    Force_Poly_One = Rom_One.Force_Polynomial;
    Stiffness_Poly_One = Rom_One.Reduced_Stiffness_Polynomial;
    
    force_equation = @(x,F) static_solution(x,F,Force_Poly_One,Stiffness_Poly_One);
    % energy_equation = @(x,lambda,F) energy_solution(x,lambda,F,energy_limit,Potential_Poly_One,Force_Poly_One,Stiffness_Poly_One);
    
    Disp_Poly_One = Rom_One.Condensed_Displacement_Polynomial;
    Disp_Poly_Two = Rom_Two.Condensed_Displacement_Polynomial;

    num_coeffs_one = Disp_Poly_One.num_element_coefficients;
    num_coeffs_two = Disp_Poly_Two.num_element_coefficients;
    input_index_two = Disp_Poly_Two.input_order;

    beta_one = Disp_Poly_One.coefficients;
    beta_two = Disp_Poly_Two.coefficients;

    beta_bar_one = beta_one*linear_mass*beta_one';
    beta_bar_two = beta_two*linear_mass*beta_two';


    Potential_Poly_Two = Rom_Two.Potential_Polynomial;

    new_sep_id = cell(1,num_validated_seps);
    new_loads = cell(1,num_validated_seps);
    new_error = cell(1,num_validated_seps);

    force_converged = zeros(1,num_validated_seps);
    disp_converged = zeros(1,num_validated_seps);

    max_sep_force_error = zeros(1,num_validated_seps);
    max_sep_disp_error = zeros(1,num_validated_seps);


    for iSep = 1:num_validated_seps %parallelise 
        estimated_end_force = scaled_force_ratios(:,iSep);
        x_0 = zeros(num_modes,1);
        
        [estimated_sep_end,~,exit_flag] = fsolve(@(x)force_equation(x,estimated_end_force),x_0,solve_opts); %#ok<ASGLU>
        
        % estimated_potential = evaluate_polynomial(Potential_Poly_One,estimated_sep_end);
        % estimated_lambda = energy_limit/estimated_potential;
        % x_0 = [estimated_sep_end;estimated_lambda];
        % [sep_end_condition,obj_value,exit_flag,output] = fsolve(@(x)energy_equation(x(1:num_modes,:),x(num_modes+1,:),estimated_end_force),x_0,solve_opts);
        
        lambda_end = find_sep_end(Potential_Poly_One,fitting_energy_limit,estimated_sep_end,estimated_end_force,1,force_equation,solve_opts);

        sep_lambda = linspace(0,lambda_end,max_sep_loadcases+1); 
        force_one = estimated_end_force.*sep_lambda(2:end);
        
        sep_r = zeros(num_modes,max_sep_loadcases);
        x_0 = zeros(num_modes,1);
        for iLoad = 1:max_sep_loadcases
            sep_r(:,iLoad) = fsolve(@(x)force_equation(x,force_one(:,iLoad)),x_0,solve_opts);
            x_0 = sep_r(:,iLoad);
        end

        
        if ismember(validated_seps(iSep),Static_Data.static_equilibrium_path_id)
            sep_distance = sqrt(sum(sep_r.^2,1));
            valid_point_spacing = false(1,max_sep_loadcases);

            sep_r_fom = Static_Data.reduced_displacement(:,Static_Data.static_equilibrium_path_id == validated_seps(iSep));
            sep_distance_fom = sqrt(sum(sep_r_fom.^2,1));

            max_distance = max(sep_distance_fom);
            min_spacing = max_distance/(1.5*max_sep_loadcases);

            

            %find closest FOM point to each ROM point
            for iPoint = 1:max_sep_loadcases
                [~,closest_point] = min(abs(sep_distance_fom-sep_distance(iPoint)));
                %check if its too close
                valid_point_spacing(iPoint) = abs(sep_distance(iPoint)-sep_distance_fom(closest_point)) > min_spacing;
            end

            sep_r = sep_r(:,valid_point_spacing);
            force_one = force_one(:,valid_point_spacing);
        end

        num_sep_points = size(sep_r,2);
        if num_sep_points == 0
            continue
        end

        potential_one = evaluate_polynomial(Potential_Poly_One,sep_r);
        potential_two = evaluate_polynomial(Potential_Poly_Two,sep_r);
        force_error = discrepency_error(potential_one,potential_two);
        
        norm_force_error = max(force_error,[],1)/max_interpolation_error(1);

        energy_one = evaluate_polynomial(Potential_Poly_One,sep_r);
        in_energy_limit = energy_one < energy_limit;
        
        if all(norm_force_error(in_energy_limit) < 1)
            force_converged(1,iSep) = 1;
        end
        max_sep_force_error(1,iSep) = max([0,norm_force_error(in_energy_limit)]);
        

        disp_prod_one = zeros(1,num_sep_points);
        disp_prod_two = zeros(1,num_sep_points);
        for iPoint = 1:num_sep_points
            r_power_products_two = ones(num_coeffs_two,1);
            for iMode = 1:num_modes
                r_power_products_two = r_power_products_two.*sep_r(iMode,iPoint).^input_index_two(:,iMode);
            end
            r_power_products_one = r_power_products_two(1:num_coeffs_one,1);


            disp_prod_one(iPoint) = r_power_products_one'*beta_bar_one*r_power_products_one;
            disp_prod_two(iPoint) = r_power_products_two'*beta_bar_two*r_power_products_two;
        end
        


        disp_error = discrepency_error(disp_prod_one,disp_prod_two);
        % [max_disp_error,max_error_index] = max(disp_error);
        norm_disp_error = disp_error/max_interpolation_error(2)^2;

        if all(norm_disp_error(in_energy_limit) < 1)
            disp_converged(1,iSep) = 1;
        end
        max_sep_disp_error(1,iSep) = max([0,norm_disp_error(in_energy_limit)]);

        interpolation_error = max(norm_force_error,norm_disp_error);
        interpolation_error(interpolation_error < 1) = 0;
        %pick worst points
        sorted_error = sort(interpolation_error,"descend");
        worst_errors = sorted_error(1:num_sep_loadcases);
        worst_errors = worst_errors(worst_errors ~= 0);
        error_index = ismember(interpolation_error,worst_errors);
        


        if any(error_index)
            num_extra_points = sum(error_index);
            new_sep_id{1,iSep} = ones(1,num_extra_points)*validated_seps(iSep);
            new_loads{1,iSep} = force_one(:,error_index);
            new_error{1,iSep} = interpolation_error(error_index);
        end

    end


    if all(force_converged == 1)
        force_log_message = "energy converged";
    else
        force_log_message = sprintf("Max energy error: %.1f" ,max(max_sep_force_error));
    end

    if all(disp_converged == 1)
        disp_log_message = "displacement converged";
    else
        disp_log_message = sprintf("Max disp error: %.1f" ,max(max_sep_disp_error));

    end

    new_sep_id = [new_sep_id{1,:}];
    new_loads = [new_loads{1,:}];
    new_error = [new_error{1,:}];

    [~,sort_index] = sort(new_error,"descend");
    
    num_extra_points = min(max_iteration_loadcases,size(new_sep_id,2));
    added_point_index = sort_index(1:num_extra_points);
    new_loads = new_loads(:,added_point_index);
    new_sep_id = new_sep_id(:,added_point_index);


    check_interpolation_time = toc(check_interpolation_start);
    interpolation_log_message = sprintf("Interpolation error check: %.1f seconds" ,check_interpolation_time);
    

    if iIteration < max_iterations + 1
        logger(force_log_message,3)
        logger(disp_log_message,3)
        logger(interpolation_log_message,3)

        if ~isempty(new_sep_id)
            if data_available
                [found_loadcases,found_r,found_theta,found_f,found_E,found_additional_data] = Static_Data_Added.contains_loadcase(new_loads);
                new_loads = new_loads(:,~found_loadcases);
            end
            [r,theta,f,E,additional_data] = Rom_One.Model.add_point(new_loads,Static_Data.additional_data_type);
            if data_available
                r(:,~found_loadcases) = r;
                r(:,found_loadcases) = found_r;

                theta(:,~found_loadcases) = theta;
                theta(:,found_loadcases) = found_theta;

                f(:,~found_loadcases) = f;
                f(:,found_loadcases) = found_f;

                E(:,~found_loadcases) = E;
                E(:,found_loadcases) = found_E;

                switch Static_Data.additional_data_type
                    case {"stiffness","perturbation"}
                        additional_data(:,:,~found_loadcases) = additional_data;
                        additional_data(:,:,found_loadcases) = found_additional_data;
                    case "none"
                end

            end
            Static_Data = Static_Data.update_data(r,theta,f,E,new_sep_id-num_original_seps,additional_data);
        end



        validation_iteration_time = toc(validation_iteration_start);
        log_message = sprintf("Validation step %i/%i completed: %i points added in %.1f seconds" ,[iIteration,max_iterations,num_extra_points,validation_iteration_time]);
        logger(log_message,2)

        if isempty(new_sep_id)
            break
        end
    else
        logger(force_log_message,2)
        logger(disp_log_message,2)
        logger(interpolation_log_message,3)
    end
end

Static_Data.validated_degree = [force_degree,disp_degree];
end

function [force,jacobian] = static_solution(x,F,Force_Poly_One,Stiffness_Poly_One)
force = Force_Poly_One.evaluate_polynomial(x) - F;
jacobian =  Stiffness_Poly_One.evaluate_polynomial(x);

end

function [sep_condition,jacobian] = energy_solution(x,lambda,F,E,Potential_Poly_One,Force_Poly_One,Stiffness_Poly_One)
potential = Potential_Poly_One.evaluate_polynomial(x) - E;
force =  Force_Poly_One.evaluate_polynomial(x) - lambda*F;
sep_condition = [potential;force];

num_modes = size(x,1);
jacobian = zeros(num_modes+1);
jacobian(1,1:num_modes) = force';
jacobian(1+(1:num_modes),1:num_modes) = Stiffness_Poly_One.evaluate_polynomial(x);
jacobian(1+(1:num_modes),1+num_modes) = -F;
end

function lambda = find_sep_end(V_Poly,V_lim,sep_end,unit_force,lambda,force_equation,solve_opts)
lambda_start = lambda;
sep_end_start = sep_end;

MAX_INC = 100;
CONVERGENCE_TOL = 0.01;
NUM_POINTS = 5;
num_modes = size(unit_force,1);


%find itial bound
sep_energy = V_Poly.evaluate_polynomial(sep_end);
energy_diff = sep_energy/V_lim;

for iInc = 1:MAX_INC 
    if abs(energy_diff - 1) < CONVERGENCE_TOL
        return
    end

    if energy_diff < 0.95
        break
    end


    exit_flag = -1;
    while exit_flag < 0
        lambda = lambda*0.95;
        [sep_end,~,exit_flag] = fsolve(@(x)force_equation(x,unit_force*lambda),zeros(num_modes,1),solve_opts);
        if lambda < 0.01
            error("Cannot find SEP solution")
        end
    end
    sep_energy = V_Poly.evaluate_polynomial(sep_end);
    energy_diff = sep_energy/V_lim;
end
if iInc == MAX_INC
    error("Couldn't find SEP end point")
end

for iInc = 1:MAX_INC
    lambda_range = linspace(lambda,lambda/energy_diff,NUM_POINTS);
    x_0 = sep_end(:,1);
    sep_end = zeros(num_modes,NUM_POINTS);
    sep_end(:,1) = x_0;
    exit_flag = zeros(1,NUM_POINTS);
    for iPoint = 2:NUM_POINTS
        [sep_end(:,iPoint),~,exit_flag(iPoint)] = fsolve(@(x)force_equation(x,unit_force*lambda_range(iPoint)),x_0,solve_opts);
        x_0 = sep_end(:,iPoint);
    end

    energy_condition = V_Poly.evaluate_polynomial(sep_end) - V_lim;
    if any(sign(energy_condition) == -1) && any(sign(energy_condition) == 1)
        break
    end

    if energy_diff < 1
        lambda = lambda*1.05;
    else
        lambda = lambda*0.95;
    end
end

if iInc == MAX_INC
    error("Couldn't find SEP end point")
end

for iInc = 1:MAX_INC

    for iPoint = 1:(NUM_POINTS-1)
        if energy_condition(iPoint)*energy_condition(iPoint + 1) < 0
            bound_indices = [iPoint,iPoint + 1];
            if iInc == 1 && any(exit_flag(bound_indices) < 0)
                warning("could not converge on bounding solutions")
                lambda = lambda_start;
                return
            end
            break
        end
    end
    lambda = interp1(energy_condition(bound_indices),lambda_range(bound_indices),0);
    x_0 = sep_end(:,bound_indices(1));
    [end_sol,~,exit_flag] = fsolve(@(x)force_equation(x,unit_force*lambda),x_0,solve_opts);
    if exit_flag < 1
        warning("solver error")
    end
    sep_energy = V_Poly.evaluate_polynomial(end_sol);
    energy_diff = sep_energy/V_lim;
    if abs(energy_diff - 1) < CONVERGENCE_TOL
        return
    end

    lambda_range = linspace(lambda_range(bound_indices(1)),lambda_range(bound_indices(2)),NUM_POINTS);
    x_0 = sep_end(:,bound_indices(1));
    x_n = sep_end(:,bound_indices(2));
    sep_end = zeros(num_modes,NUM_POINTS);
    sep_end(:,1) = x_0;
    sep_end(:,end) = x_n;
    for iPoint = 2:(NUM_POINTS-1)
        [sep_end(:,iPoint),~,exit_flag] = fsolve(@(x)force_equation(x,unit_force*lambda_range(iPoint)),x_0,solve_opts);
        if exit_flag < 1
            warning("solver error")
        end
        x_0 = sep_end(:,iPoint);
    end

    energy_condition = V_Poly.evaluate_polynomial(sep_end) - V_lim;
end

if iInc == MAX_INC
    warning("Couldn't find SEP end point: V_end/V_lim = " + energy_diff)
end

end