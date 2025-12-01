function Static_Data = verify_validation_polynomials(Static_Data)
MAXIMUM_DEGREE = 12;
load_data = 1;
save_data = 1;

Verification_Opts = Static_Data.Verification_Options;
max_interpolation_error = Verification_Opts.maximum_interpolation_error(3:4);

Model = Static_Data.Model;
max_sep_points = Model.Static_Options.maximum_loadcases;

initial_degree = Static_Data.verified_degree;
% initial_degree = Static_Data.verified_degree - 1;



num_r_modes = length(Model.reduced_modes);

adjacent_sep_ratios = get_adjacent_sep_ratios(Static_Data.unit_sep_ratios);

unit_force_ratios = [Static_Data.unit_sep_ratios,adjacent_sep_ratios];
[~,~,reorder_index] = uniquetol(unit_force_ratios', "ByRows",1);
removal_indicies = false(1,size(unit_force_ratios,2));
for iIndex = 1:max(reorder_index)
    matching_indices = reorder_index == iIndex;
    if nnz(matching_indices) > 1
        repeated_index = find(matching_indices);
        repeated_index(1) = [];
        removal_indicies(repeated_index) = 1;
    end
end
unit_force_ratios(:,removal_indicies) = [];

scaled_force_ratios = scale_sep_ratios(unit_force_ratios,Static_Data.Model.calibrated_forces);
num_verified_seps = size(unit_force_ratios,2);


num_dataset_points = size(Static_Data,2);
max_degree = get_max_poly_degree("displacement",num_r_modes,num_dataset_points,MAXIMUM_DEGREE);
% max_degree = max_degree - 1;
%-------------
stiffness_converged = zeros(1,num_verified_seps);
disp_grad_converged = zeros(1,num_verified_seps);

degree = initial_degree;
Static_Data.Dynamic_Validation_Data.degree = degree;
Rom_One = Reduced_System(Static_Data,"id",2);

num_degree_pairs = min(max_degree+2 - degree(1),max_degree+2 - degree(2))/2;
maximum_stiffness_pair_errors = zeros(1,num_degree_pairs);
maximum_disp_grad_pair_errors = zeros(1,num_degree_pairs);
max_pair_error = inf;

num_validation_modes = size(Static_Data.low_frequency_stiffness,1);
validation_points = add_sep_ratios(num_validation_modes,2);

fitting_energy_limit = Model.fitting_energy_limit;
energy_limit = Model.energy_limit;


verification_data = cell(1,num_degree_pairs);
error_calculation_failed = 0;
error_time_start = tic;

data_path = Static_Data.get_data_path + "validation\sep_rom_data.mat";
if load_data
    data_available = isfile(data_path);
    if data_available
        load(data_path,"sep_rom_data")
        if size(sep_rom_data,2) ~= num_verified_seps
            data_available = false;
        end
    end
else
    empty_cell = cell(1,num_verified_seps);
    sep_rom_data = struct("disp",empty_cell,"lambda",empty_cell);
    data_available = false;
end

for iDegree_pair = 1:num_degree_pairs
    
    maximum_stiffness_pair_error = 0;
    maximum_disp_grad_pair_error = 0;

    degree_two = degree + 2;



    log_message = sprintf("Comparing %s and %s degree stiffness polynomials" , ...
        ordinal_suffix(degree(1)),ordinal_suffix(degree_two(1)));
    logger(log_message,4)
    log_message = sprintf("Comparing %s and %s degree displacement gradient polynomials" , ...
        ordinal_suffix(degree(2)),ordinal_suffix(degree_two(2)));
    logger(log_message,4)
    
    Static_Data.Dynamic_Validation_Data.degree = degree_two;
    Rom_Two = Reduced_System(Static_Data,"id",3);


    Disp_Error_Inputs.Beta_Bar_Data_One = Rom_One.get_h_beta_bar(Rom_One.Low_Frequency_Coupling_Gradient_Polynomial.coefficients,Rom_One.Physical_Displacement_Polynomial.coefficients);
    Disp_Error_Inputs.Beta_Bar_Data_Two = Rom_Two.get_h_beta_bar(Rom_Two.Low_Frequency_Coupling_Gradient_Polynomial.coefficients,Rom_Two.Physical_Displacement_Polynomial.coefficients);
    Disp_Error_Inputs.input_order = Rom_Two.get_max_input_order;
    new_error = cell(1,num_verified_seps);


    log_message = sprintf("Checking %d SEPs..." ,num_verified_seps);
    logger(log_message,4)

    error_calculation_failed = zeros(1,num_verified_seps);
    load("data\plot_level.mat","plotting_level")
    if plotting_level >= 4
        num_jobs = 0;
        warning("parallelisation disabled for verification plotting")
    else
        num_jobs = get_current_parallel_jobs;
    end




    Disp_Error_Inputs_Const = parallel.pool.Constant(Disp_Error_Inputs);
    Rom_One_Const =  parallel.pool.Constant(Rom_One);
    Rom_Two_Const = parallel.pool.Constant(Rom_Two);
    parfor (iSep = 1:num_verified_seps,num_jobs)
    % for iSep = 1:num_verified_seps
         
        force_ratio = scaled_force_ratios(:,iSep);
        if load_data && data_available
            disp_sep = sep_rom_data(iSep).disp;
            lambda_sep = sep_rom_data(iSep).lambda;
        else
            [disp_sep,lambda_sep] = find_sep_rom(Rom_One_Const.Value,force_ratio,3*max_sep_points);
            if save_data
                sep_rom_data(iSep).disp = disp_sep;
                sep_rom_data(iSep).lambda = lambda_sep;
            end
        end
        
        if isempty(lambda_sep)
            error_calculation_failed(iSep) = 1;
        else
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

            stiffness_error = get_stiffness_error(tested_disp,validation_points,Rom_One_Const.Value,Rom_Two_Const.Value);
            disp_grad_error = get_disp_gradient_error(tested_disp,validation_points,Rom_One_Const.Value,Rom_Two_Const.Value,force_ratio,Disp_Error_Inputs_Const.Value);

            %---------------


            norm_stiffness_error = stiffness_error/max_interpolation_error(1);
            norm_disp_grad_error = disp_grad_error/max_interpolation_error(2);

            potential_one = Rom_One_Const.Value.Potential_Polynomial.evaluate_polynomial(tested_disp);
            in_energy_limit = potential_one <= energy_limit;

            % Plot_Data = struct([]);
            % Plot_Data(1).displacement = disp_sep;
            % Plot_Data(1).lambda = lambda_sep;
            % Plot_Data(1).Rom = Rom_One_Const.Value;
            % Plot_Data(1).force_ratio = force_ratio;
            % Plot_Data(1).tested_index = tested_sep_index;
            % Plot_Data(2).Rom = Rom_Two_Const.Value;
            % 
            % Plot_Error = struct("force_error",[],"disp_error",[]);
            % Plot_Error.force_error = norm_stiffness_error;
            % Plot_Error.disp_error = norm_disp_grad_error;
            % 
            % if max(norm_stiffness_error) > 10
            %     figure
            %     tiledlayout("flow")
            %     nexttile
            %     plot(lambda_sep(tested_sep_index),norm_stiffness_error)
            %     nexttile
            %     plot(lambda_sep(tested_sep_index),norm_disp_grad_error)
            %     1
            % end
            % sep_error_plot(Plot_Data,Static_Data,Plot_Error,Disp_Error_Inputs_Const.Value)

            if all(norm_stiffness_error(in_energy_limit) < 1)
                stiffness_converged(1,iSep) = 1;
            end
            if all(norm_disp_grad_error(in_energy_limit) < 1)
                disp_grad_converged(1,iSep) = 1;
            end


            maximum_stiffness_pair_error = max(maximum_stiffness_pair_error,max(norm_stiffness_error));
            maximum_disp_grad_pair_error = max(maximum_disp_grad_pair_error,max(norm_disp_grad_error));

            interpolation_error = max(norm_stiffness_error,norm_disp_grad_error);
            interpolation_error(interpolation_error < 1) = 0;



            %pick worst points
            [sorted_error,sorted_error_index] = sort(interpolation_error,"descend");
            % num_error_points = size(sorted_error,2);

            worst_errors = sorted_error(1);
            worst_errors_index = sorted_error_index(1);

            % for auto, add the worst point and
            %   - for 3 add the error = 1 points


            nonzero_error_index = worst_errors ~= 0;
            worst_errors = worst_errors(nonzero_error_index);
            error_index = worst_errors_index(nonzero_error_index);

     
        end
    end
    clear("Disp_Error_Inputs_Const","Rom_One_Const","Rom_Two_Const")


    Extra_Point_Data.new_error = new_error;
    Extra_Point_Data.degree = degree;
    verification_data{1,iDegree_pair} = Extra_Point_Data;


    maximum_stiffness_pair_errors(iDegree_pair) = maximum_stiffness_pair_error;
    maximum_disp_grad_pair_errors(iDegree_pair) = maximum_disp_grad_pair_error;

    log_message = sprintf("Max stiffness error: %.2f and max disp grad error: %.2f" ,maximum_stiffness_pair_error,maximum_disp_grad_pair_error);
    logger(log_message,4)

    %if the error starts going up terminate early
    % last_max_pair_error = max_pair_error;
    % max_pair_error = max(maximum_disp_pair_errors(iDegree_pair), maximum_force_pair_errors(iDegree_pair));
    % if max_pair_error > last_max_pair_error
    %     maximum_disp_pair_errors(:,(iDegree_pair+1):end) = [];
    %     maximum_force_pair_errors(:,(iDegree_pair+1):end) = [];
    %     break
    % end


    if ~all(stiffness_converged) || ~all(disp_grad_converged)
        degree(1) = degree_two(1);
    end

    if ~all(disp_grad_converged)
        degree(2) = degree_two(2);
    end

    if all(stiffness_converged) && all(disp_grad_converged)
        maximum_disp_grad_pair_errors(:,(iDegree_pair+1):end) = [];
        maximum_stiffness_pair_errors(:,(iDegree_pair+1):end) = [];
        break
    end
    
    Static_Data.Dynamic_Validation_Data.degree = degree;
    Rom_One = Reduced_System(Static_Data,"id",2);
end
if save_data && ~data_available
    dir_path = Static_Data.get_data_path + "validation";
    if ~isfolder(dir_path)
        mkdir(dir_path)
    end
    save(data_path,"sep_rom_data")
end
error_time = toc(error_time_start);
%------------------------------------------------------------------
if any(error_calculation_failed)
    log_message = sprintf("Error calculation failed. Could not follow SEPs");
    logger(log_message,1)
end
[~,stiffness_degree_index] = min(maximum_stiffness_pair_errors);
[~,disp_grad_degree_index] = min(maximum_disp_grad_pair_errors);

degree(1) = verification_data{1,stiffness_degree_index}.degree(1);
degree(2) = verification_data{1,disp_grad_degree_index}.degree(2);

Static_Data.Dynamic_Validation_Data.degree = degree;
log_message = sprintf("Max stiffness error: %.2f and max displacement gradient error: %.2f found in %.1f seconds" ,maximum_stiffness_pair_errors(stiffness_degree_index),maximum_disp_grad_pair_errors(disp_grad_degree_index),error_time);
    logger(log_message,3)
end