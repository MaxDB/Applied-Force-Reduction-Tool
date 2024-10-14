function Static_Data = sep_verification(Static_Data)
MAXIMUM_DEGREE = 13; %maximum degree that can be fit
INITIAL_FORCE_DEGREE = 3;
INITIAL_DISPLACEMENT_DEGREE = 3;

%-------------------------
Verification_Opts = Static_Data.Verification_Options;
max_interpolation_error = Verification_Opts.maximum_interpolation_error;
max_iteration_loadcases_setting = Verification_Opts.max_added_points;
max_iterations = Verification_Opts.maximum_iterations;

num_added_points = Verification_Opts.num_added_points;

Model = Static_Data.Model;
max_sep_points = Model.Static_Options.maximum_loadcases;

fitting_energy_limit = Model.fitting_energy_limit;
energy_limit = Model.energy_limit;

num_r_modes = length(Model.reduced_modes);


%-------------------------
found_force_ratios = Static_Data.unit_sep_ratios;
num_original_seps = size(found_force_ratios,2);



unit_force_ratios = found_force_ratios;

scaffold_points = Static_Data.scaffold_points;
data_available = any(scaffold_points == 0);
if data_available
    added_data = ~scaffold_points;
    [Static_Data,Static_Data_Added] = Static_Data.remove_loadcases(added_data);
end

force_degree = INITIAL_FORCE_DEGREE;
disp_degree = INITIAL_DISPLACEMENT_DEGREE;
sep_density = 1;
for iIteration = 1:(max_iterations+1)
    sep_density = sep_density + 2;
    new_unit_force_ratios = add_sep_ratios(num_r_modes,sep_density,found_force_ratios);
    unit_force_ratios = [unit_force_ratios,new_unit_force_ratios]; %#ok<AGROW>

    scaled_force_ratios = scale_sep_ratios(unit_force_ratios,Static_Data.Model.calibrated_forces);
    found_force_ratios = unit_force_ratios;
    num_verified_seps = size(unit_force_ratios,2);
    
    validated_seps = (1:num_verified_seps);

    V = Static_Data.potential_energy;
    num_dataset_points = size(V,2);

    max_force_degree = get_max_poly_degree("force",num_r_modes,num_dataset_points,MAXIMUM_DEGREE);
    max_disp_degree = get_max_poly_degree("displacement",num_r_modes,num_dataset_points,MAXIMUM_DEGREE);
    max_degree = max(max_force_degree,max_disp_degree);
    next_max_degree = min(MAXIMUM_DEGREE,max_degree + 2);
    data_size = size(Static_Data);
    max_iteration_loadcases = get_max_added_points(next_max_degree,data_size,...
        max_iteration_loadcases_setting);

    if num_added_points == 0
        continue
    end
    
    force_converged = zeros(1,num_verified_seps);
    disp_converged = zeros(1,num_verified_seps);

    Rom_One = Reduced_System(Static_Data,[force_degree,disp_degree]);
    
    num_degree_pairs = max(max_force_degree - force_degree,max_disp_degree - disp_degree)/2;
    maximum_force_pair_errors = zeros(1,num_degree_pairs);
    maximum_disp_pair_errors = zeros(1,num_degree_pairs);
    max_pair_error = inf;

    verification_data = cell(1,num_degree_pairs);
    error_time_start = tic;
    for iDegree_pair = 1:num_degree_pairs
        maximum_force_pair_error = 0;
        maximum_disp_pair_error = 0;

        force_degree_two = force_degree + 2;
        disp_degree_two = disp_degree + 2;


        Rom_Two = Reduced_System(Static_Data,[force_degree_two,disp_degree_two]);

        Disp_Error_Inputs.beta_bar_one = Rom_One.get_beta_bar(Rom_One.Physical_Displacement_Polynomial);
        Disp_Error_Inputs.beta_bar_two = Rom_Two.get_beta_bar(Rom_Two.Physical_Displacement_Polynomial);
        Disp_Error_Inputs.input_order = Rom_Two.get_max_input_order;

        Disp_Error_Inputs.Disp_Diff_Data_One = Rom_One.Physical_Displacement_Polynomial.get_diff_data(1);
        Disp_Error_Inputs.Disp_Diff_Data_Two = Rom_Two.Physical_Displacement_Polynomial.get_diff_data(1);

        

        max_sep_force_error = zeros(1,num_verified_seps);
        max_sep_disp_error = zeros(1,num_verified_seps);

        new_sep_id = cell(1,num_verified_seps);
        new_loads = cell(1,num_verified_seps);
        new_error = cell(1,num_verified_seps);
        
        for iSep = 1:num_verified_seps
            validation_iteration_start = tic;

            force_ratio = scaled_force_ratios(:,iSep);
            [disp_sep,lambda_sep] = find_sep_rom(Rom_One,force_ratio,3*max_sep_points);
            lambda_step = lambda_sep(end)/max_sep_points;

            tested_sep_index = zeros(1,max_sep_points);
            norm_lambda_sep = round(lambda_sep/lambda_step);
            next_sep_index = 1;
            for iPoint = 1:max_sep_points
                next_sep_index = find(norm_lambda_sep(1,next_sep_index:end) >= iPoint,1) + next_sep_index - 1;
                tested_sep_index(iPoint) = next_sep_index;
            end

            force_sep = lambda_sep.*force_ratio;
            %----------------

            tested_disp = disp_sep(:,tested_sep_index);
            tested_force = force_sep(:,tested_sep_index);
            
            force_error = get_force_error(tested_disp,Rom_One,Rom_Two);
            disp_error = get_disp_error(tested_disp,Rom_One,Rom_Two,force_ratio,Disp_Error_Inputs);

            %---------------
            Potential_One = Rom_One.Potential_Polynomial.evaluate_polynomial(tested_disp);
            in_energy_limit = Potential_One <= energy_limit;

            norm_force_error = force_error/max_interpolation_error(1);
            norm_disp_error = disp_error/max_interpolation_error(2);
            
            if all(norm_force_error(in_energy_limit) < 1)
                force_converged(1,iSep) = 1;
            end
            if all(norm_disp_error(in_energy_limit) < 1)
                disp_converged(1,iSep) = 1;
            end
            max_sep_force_error(1,iSep) = max(norm_force_error);
            max_sep_disp_error(1,iSep) = max(norm_disp_error);
            
            maximum_force_pair_error = max(maximum_force_pair_error,max(norm_force_error));
            maximum_disp_pair_error = max(maximum_disp_pair_error,max(norm_disp_error));

            interpolation_error = max(norm_force_error,norm_disp_error);
            interpolation_error(interpolation_error < 1) = 0;

            

            %pick worst points
            [sorted_error,sorted_error_index] = sort(interpolation_error,"descend");
            num_error_points = size(sorted_error,2);

            worst_errors = sorted_error(1:min(num_error_points,num_added_points));
            worst_errors_index = sorted_error_index(1:min(num_error_points,num_added_points));

            nonzero_error_index = worst_errors ~= 0;
            worst_errors = worst_errors(nonzero_error_index);
            error_index = worst_errors_index(nonzero_error_index);

            if ~isempty(error_index)
                num_extra_points = size(error_index,2);
                new_sep_id{1,iSep} = ones(1,num_extra_points)*validated_seps(iSep);
                new_loads{1,iSep} = tested_force(:,error_index);
                new_error{1,iSep} = worst_errors;
            end
        end


        Extra_Point_Data.new_sep_id = new_sep_id;
        Extra_Point_Data.new_loads = new_loads;
        Extra_Point_Data.new_error = new_error;
        Extra_Point_Data.degree = [force_degree,disp_degree];
        verification_data{1,iDegree_pair} = Extra_Point_Data;


        maximum_force_pair_errors(iDegree_pair) = maximum_force_pair_error;
        maximum_disp_pair_errors(iDegree_pair) = maximum_disp_pair_error;

        %if the error starts going up terminate early
        last_max_pair_error = max_pair_error;
        max_pair_error = max(maximum_disp_pair_errors(iDegree_pair), maximum_force_pair_errors(iDegree_pair));
        if max_pair_error > last_max_pair_error
            maximum_disp_pair_errors(:,(iDegree_pair+1):end) = [];
            maximum_force_pair_errors(:,(iDegree_pair+1):end) = [];
            break
        end


        if ~all(force_converged)
            force_degree = force_degree_two;
        end

        if ~all(disp_converged)
            disp_degree = disp_degree_two;
        end

        if all(force_converged) && all(disp_converged)
            maximum_disp_pair_errors(:,(iDegree_pair+1):end) = [];
            maximum_force_pair_errors(:,(iDegree_pair+1):end) = [];
            break
        end
        
        Rom_One = Reduced_System(Static_Data,[force_degree,disp_degree]);
        % if all(force_converged) || all(disp_converged)
        %     
        % else
        %     Rom_One = Rom_Two;
        % end
    end
    error_time = toc(error_time_start);
            %------------------------------------------------------------------
    [~,degree_index] = min(max(maximum_disp_pair_errors,maximum_force_pair_errors)); 
    Extra_Point_Data = verification_data{1,degree_index};
    force_degree = Extra_Point_Data.degree(1);
    disp_degree = Extra_Point_Data.degree(2);

    log_message = sprintf("Max force error: %.2f and max disp error: %.2f found in %.1f seconds" ,maximum_force_pair_errors(degree_index),maximum_disp_pair_errors(degree_index),error_time);
    logger(log_message,3)

    new_sep_id = [Extra_Point_Data.new_sep_id{1,:}];
    new_loads = [Extra_Point_Data.new_loads{1,:}];
    new_error = [Extra_Point_Data.new_error{1,:}];

    [~,sort_index] = sort(new_error,"descend");

    num_extra_points = min(max_iteration_loadcases,size(new_sep_id,2));
    added_point_index = sort_index(1:num_extra_points);
    new_loads = new_loads(:,added_point_index);
    new_sep_id = new_sep_id(:,added_point_index);

    if iIteration > max_iterations
        break
    end

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

        removal_index = E > fitting_energy_limit;
        if any(removal_index)
            r(:,removal_index) = [];
            theta(:,removal_index) = [];
            f(:,removal_index) = [];
            E(:,removal_index) = [];
            new_sep_id(:,removal_index) = [];
            switch Static_Data.additional_data_type
                case {"stiffness","perturbation"}
                    additional_data(:,:,removal_index) = [];
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


end

Static_Data.verified_degree = [force_degree,disp_degree];
end

function max_iteration_loadcases = get_max_added_points(degree,data_size,max_loadcase_setting)
num_r_modes = data_size(1);
num_points = data_size(2);
min_points = 0;
for iDegree = 2:degree
    min_points = min_points + nchoosek(iDegree + num_r_modes - 1,num_r_modes-1);
end

point_diff = min_points - num_points;
max_iteration_loadcases = max(point_diff,max_loadcase_setting);
end