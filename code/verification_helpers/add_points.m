function Static_Data = add_points(seps,force_ratios,Static_Data,rom_one,rom_two,error_modes,max_interpolation_error)
V = Static_Data.potential_energy;
energy_span = V <= rom_one.Model.energy_limit;
max_sep_points = rom_one.Model.Static_Options.maximum_loadcases;

r_disp = Static_Data.reduced_displacement(:,energy_span);
f_r = Static_Data.restoring_force(:,energy_span);
sep_id = Static_Data.static_equilibrium_path_id(:,energy_span);

num_r_modes = size(r_disp,1);
num_seps = length(seps);
new_sep_id = [];
new_loads = zeros(num_r_modes,0);
num_restarted_seps = 0;
for iSep = 1:num_seps
    sep = seps(iSep);
    error_mode = error_modes(iSep);

    sep_span = sep_id == sep;
    num_sep_points = sum(sep_span);

    if num_sep_points >= max_sep_points
        continue
    end

    f_sep = f_r(:,sep_span);
    [~,sep_order] = sort(max(abs(f_sep),[],1),"ascend");
    f_sep = f_sep(:,sep_order);

    [disp_sep_one,force_sep_one] = find_sep_rom(rom_one,force_ratios(:,iSep));

    f_sep_mid = f_sep - diff([zeros(num_r_modes,1),f_sep],1,2)/2;
    sep_force_index = find(f_sep_mid(:,1) ~= 0,1);
    zero_force_index = (f_sep_mid(:,1) == 0);

    r_sep_mid = zeros(size(f_sep_mid));
    for iMode = 1:num_r_modes
        r_sep_mid(iMode,:) = interp1(force_sep_one(sep_force_index,:),disp_sep_one(iMode,:),f_sep_mid(sep_force_index,:));
    end

    f_one = rom_one.Force_Polynomial.evaluate_polynomial(r_sep_mid);
    f_two = rom_two.Force_Polynomial.evaluate_polynomial(r_sep_mid);
    force_error = discrepency_error(f_one,f_two);
    force_error(zero_force_index,:) = 0;
    max_force_error = max(force_error,[],1);

    theta_one = rom_one.Condensed_Displacement_Polynomial.evaluate_polynomial(r_sep_mid,error_mode);
    theta_two = rom_two.Condensed_Displacement_Polynomial.evaluate_polynomial(r_sep_mid,error_mode);
    theta_error = discrepency_error(theta_one,theta_two);

    add_sep_points = max_force_error > max_interpolation_error(1) | theta_error > max_interpolation_error(2);

    num_extra_points = sum(add_sep_points);
    if num_sep_points + num_extra_points > max_sep_points
        theta_error_points = theta_error/max_interpolation_error(2);
        theta_error_points(~add_sep_points) = 0;
        
        force_error_points = max_force_error/max_interpolation_error(1);
        force_error_points(~add_sep_points) = 0;

        error_points = max(theta_error_points,force_error_points);
        error_points(error_points == 0) = inf;
        %only select greatest errors within point limit
        num_removed_extra_points = num_sep_points + num_extra_points - max_sep_points;
        for iPoint = 1:num_removed_extra_points
            [~,min_error_index] = min(error_points);  %prioritise force
            add_sep_points(min_error_index) = false;
            error_points(min_error_index) = inf;
        end
    end
    num_extra_points = sum(add_sep_points);

    if any(add_sep_points)
        num_restarted_seps = num_restarted_seps + 1;
        
        new_sep_id = [new_sep_id,ones(1,num_extra_points)*sep]; %#ok<AGROW>
        new_loads = [new_loads,f_sep_mid(:,add_sep_points)]; %#ok<AGROW>
    end
end

if ~isempty(new_sep_id)
    num_original_seps = size(Static_Data.unit_sep_ratios,2);
    log_message = sprintf("%i loadcases over %i/%i SEPs added" ,[length(new_sep_id),num_restarted_seps,num_original_seps]);
    logger(log_message,3)
    [r,theta,f,E,additional_data] = rom_one.Model.add_point(new_loads,Static_Data.additional_data_type);
    Static_Data = Static_Data.update_data(r,theta,f,E,new_sep_id-num_original_seps,additional_data);
end
end