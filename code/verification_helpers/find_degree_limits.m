function Min_Degree_Data = find_degree_limits(Model,reduced_displacement,displacement,force,energy,sep_id,seps)
MIN_POINTS = 4;

mode = seps(end)/2;
Model = Model.get_modal_subset(mode);

Calibration_Opts = Model.Calibration_Options;

sep_one_index = sep_id == seps(1);
sep_two_index = sep_id == seps(2);
sep_index = sep_one_index | sep_two_index;


if nnz(sep_one_index) >= nnz(sep_two_index)
    main_sep_index = sep_one_index;
else
    main_sep_index = sep_two_index;
end

energy_main_sep = energy(main_sep_index);
num_points = size(energy_main_sep,2);


Verification_Opts.maximum_iterations = 0;
Verification_Opts.maximum_interpolation_error = [1e-4,1e-4];
Verification_Opts.verification_algorithm = "sep_grow";

energy_mode = energy(sep_index);
r_mode = reduced_displacement(mode,sep_index);
displacement_mode = displacement(:,sep_index);
force_mode = force(mode,sep_index);
sep_id_mode = sep_id(sep_index);

min_degree = nan(num_points,2);

for iPoint = 2:num_points
    energy_lim = energy_main_sep(iPoint-1);
    if energy_lim > Model.energy_limit
        break
    end
    fitting_energy_lim = energy_main_sep(iPoint);

    
    in_limit = energy_mode <= fitting_energy_lim;
    num_sub_points = nnz(in_limit);
    if num_sub_points < MIN_POINTS
        continue
    end
    

    New_Model = Model;
    New_Model.energy_limit = energy_lim;
    New_Model.fitting_energy_limit = fitting_energy_lim;
    

    Dataset.reduced_displacement = r_mode(in_limit);
    Dataset.physical_displacement = displacement_mode(:,in_limit);
    Dataset.restoring_force = force_mode(in_limit);
    Dataset.potential_energy = energy_mode(in_limit);
    Dataset.static_equilibrium_path_id = sep_id_mode(in_limit);
    Static_Data = Static_Dataset(New_Model,Verification_Opts,Dataset);

    Static_Data = Static_Data.verify_dataset;
    min_degree(iPoint,:) = Static_Data.verified_degree;
end

not_found = isnan(min_degree(:,1));
min_system_degree = min_degree(~not_found,:);
energy_lim = energy_main_sep(find(~not_found)-1);

max_force_degree = max(min_system_degree(:,1));
max_disp_degree = max(min_system_degree(:,2));

force_degrees = 3:2:max_force_degree;
disp_degrees = 3:2:max_disp_degree;

force_degree_indices = get_last_indices(min_system_degree(:,1),force_degrees);
disp_degree_indices = get_last_indices(min_system_degree(:,2),disp_degrees);

Min_Degree_Data.force_degrees = force_degrees;
Min_Degree_Data.disp_degrees = disp_degrees;
Min_Degree_Data.force_degree_energy = energy_lim(force_degree_indices);
Min_Degree_Data.disp_degree_energy = energy_lim(disp_degree_indices);
end

function indices = get_last_indices(list,values)
num_values = length(values);
indices = zeros(num_values,1);
for iValue = 1:num_values
    indices(iValue) = find(list == values(iValue),1,"last");
end
end