function Static_Data = sep_validation(Static_Data)
MAXIMUM_DEGREE = [21,15,13,11];
%CAN BE MADE A LOT MORE FFICIENT

Validation_Opts = Static_Data.Validation_Options;
max_interpolation_error = Validation_Opts.maximum_interpolation_error;
max_fitting_error = Validation_Opts.maximum_fitting_error;
min_coupling_rating = Validation_Opts.minimum_coupling_rating;

validated_dofs = Static_Data.rank_coupling(min_coupling_rating);

initial_degree = Validation_Opts.minimum_degree;
max_iterations = Validation_Opts.maximum_iterations;

num_modes = length(Static_Data.Model.reduced_modes);
max_degree = MAXIMUM_DEGREE(num_modes);


force_degree = initial_degree;
disp_degree = initial_degree;
for iIteration = 1:(max_iterations+1)
    validation_iteration_start = tic;
    V = Static_Data.potential_energy;
    validated_points = V < Static_Data.Model.energy_limit;

    r = Static_Data.reduced_displacement(:,validated_points);
    f = Static_Data.restoring_force(:,validated_points);
    theta = Static_Data.condensed_displacement(:,validated_points);

    %find minimum degree
    minimum_force_degree_start = tic;
    while force_degree <= max_degree
        rom = Reduced_System(Static_Data,[force_degree,1],"full");
        f_rom = rom.Force_Polynomial.evaluate_polynomial(r);
        force_error = coeff_of_determination(f,f_rom,rom.Force_Polynomial.num_independent_element_coefficients);
        max_force_error = max(force_error);

        if max_force_error < max_fitting_error
            break
        end

        force_degree = force_degree + 2;
    end
    if force_degree > max_degree
        error("Force polynomial requires too high a degree")
    end
    minimum_force_degree_time = toc(minimum_force_degree_start);
    logger("---",3)
    log_message = sprintf("Minimum force degree is %i: %.1f seconds" ,[force_degree,minimum_force_degree_time]);
    logger(log_message,3)

    minimum_disp_degree_time_start = tic;
    while disp_degree <= max_degree
        rom = Reduced_System(Static_Data,[1,disp_degree],"full");
        theta_rom = rom.Condensed_Displacement_Polynomial.evaluate_polynomial(r);
        disp_error = coeff_of_determination(theta(validated_dofs,:),theta_rom(validated_dofs,:),rom.Condensed_Displacement_Polynomial.num_independent_element_coefficients);

        [max_disp_error,max_disp_error_index] = max(disp_error);

        if max_disp_error < max_fitting_error
            break
        end

        disp_degree = disp_degree + 2;
    end
    if disp_degree > max_degree
        warning("Quasi-static coupling polynomial requires too high a degree. Worst validated index: " + max_disp_error_index)
    end

    minimum_disp_degree_time = toc(minimum_disp_degree_time_start);
    log_message = sprintf("Minimum disp degree is  %i: %.1f seconds" ,[disp_degree,minimum_disp_degree_time]);
    logger(log_message,3)

    degree = [force_degree,disp_degree];
    validation_plot(1,Static_Data,rom,validated_dofs(max_disp_error_index))

    %---------------------------------------------------------------------%
    %interpolation withinin SEPs


    rom_one = Reduced_System(Static_Data,degree,"full");
    rom_two = Reduced_System(Static_Data,degree + 2,"full");
    found_force_ratios = Static_Data.unit_sep_ratios;

    scaled_found_force_ratios = scale_sep_ratios(found_force_ratios,Static_Data.Model.calibrated_forces);

    [found_sep_force_error,found_sep_disp_error,~,error_modes] = check_seps_rom(rom_one,rom_two,scaled_found_force_ratios,validated_dofs); %#ok<ASGLU>
    found_sep_force_error(isinf(found_sep_force_error)) = 0;
    add_sep_points = found_sep_force_error > max_interpolation_error(1) | found_sep_disp_error > max_interpolation_error(2);
    if any(add_sep_points)
        Static_Data = add_points(find(add_sep_points),scaled_found_force_ratios(:,add_sep_points),Static_Data,rom_one,rom_two,error_modes(add_sep_points),max_interpolation_error);
        rom_one = Reduced_System(Static_Data,degree,"full");
        rom_two = Reduced_System(Static_Data,degree + 2,"full");
    end




    if iIteration ~= (max_iterations+1)
        %---------------------------------------------------------------------%
        %interpolation between SEPs

        unit_force_ratios = add_sep_ratios(num_modes,1+2*iIteration,found_force_ratios);
        scaled_force_ratios = scale_sep_ratios(unit_force_ratios,Static_Data.Model.calibrated_forces);

        [sep_force_error,sep_disp_error,sep_limit,error_modes] = check_seps_rom(rom_one,rom_two,scaled_force_ratios,validated_dofs); %#ok<ASGLU>



        add_sep = sep_force_error > max_interpolation_error(1) | sep_disp_error > max_interpolation_error(2);
        if all(~add_sep)
            break
        end

        predicted_force_ratios = sep_limit(:,add_sep)*Static_Data.Model.Calibration_Options.force_overcalibration;

        [r,theta,f,E,sep_id,additional_data] = Static_Data.Model.add_sep(predicted_force_ratios,Static_Data.additional_data_type,1);
        found_force_ratios = unit_force_ratios(:,add_sep);
        Static_Data = Static_Data.update_data(r,theta,f,E,sep_id,additional_data,found_force_ratios);
    end
    validation_iteration_time = toc(validation_iteration_start);
    log_message = sprintf("Validation step %i/%i completed: %.1f seconds" ,[iIteration,max_iterations+1,validation_iteration_time]);
    logger(log_message,2)
end
Static_Data.validated_degree = degree;
end