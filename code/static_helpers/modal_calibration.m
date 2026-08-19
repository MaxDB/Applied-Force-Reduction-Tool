function Model = modal_calibration(Model,modes)
Calibration_Opts = Model.Calibration_Options;

%check if already calibrated
GEOMETRY_PATH = "geometry\" + Model.system_name + "\";
if isfile(GEOMETRY_PATH + "force_calibration.mat")
    load(GEOMETRY_PATH + "force_calibration.mat","Force_Calibration");
end

if isfile(GEOMETRY_PATH + "force_calibration.mat") && isequal(Force_Calibration.Parameters,Model.Parameters)
    calibrated_energy = Force_Calibration.energy_limit;
    if ismember(Model.energy_limit,calibrated_energy)
        calibration_id = find(calibrated_energy == Model.energy_limit);
        calibrated_modes = Force_Calibration.calibrated_modes{1,calibration_id};
        if any(calibrated_modes > 1000)
            nc_mode_id = calibrated_modes>1000;
            calibrated_modes(nc_mode_id) = [];
            Force_Calibration.calibrated_modes{1,calibration_id} = calibrated_modes;
            f_lim = Force_Calibration.force_limit{1,calibration_id};
            f_lim(nc_mode_id,:) = [];
            Force_Calibration.force_limit{1,calibration_id} = f_lim;
            
        end
        uncalibrated_modes = setdiff(modes,calibrated_modes);
    else
        calibration_id = length(calibrated_energy) + 1;
        Force_Calibration.energy_limit(calibration_id) = Model.energy_limit;
        uncalibrated_modes = modes;
        calibrated_modes = [];
    end
else
    uncalibrated_modes = modes;
    calibrated_modes = [];
    Force_Calibration.energy_limit = Model.energy_limit;
    Force_Calibration.force_limit = {};
    Force_Calibration.calibrated_modes = {};
    calibration_id = 1;
end

%calibrate force scale factors


num_calibrated_modes = length(calibrated_modes);
num_uncalibrated_modes = length(uncalibrated_modes);


modes = Model.reduced_modes;
num_modes = length(modes);
eigenvalues = Model.reduced_eigenvalues;

num_matching_calibrated_modes = length(intersect(modes,calibrated_modes));

log_message = sprintf("%u/%u modes precalibrated",[num_matching_calibrated_modes,num_modes]);
logger(log_message,3)

initial_force_ratio = zeros(num_modes,num_uncalibrated_modes*2);
limit_type = repmat("energy",num_uncalibrated_modes,1);
for iMode = 1:num_uncalibrated_modes
    mode = uncalibrated_modes(iMode);
    mode_index = mode == modes;

    force_ratio = zeros(num_modes,2);
    force_ratio(mode_index,:) = [1,-1];

    %start with linear approximation

    force_scale_factor = Calibration_Opts.calibration_scale_factor*sqrt(2*eigenvalues(mode_index)*Model.fitting_energy_limit);
    if mode > 1000 && ~isempty(Model.nc_amplitude_limit)
        force_amp_scale_factor = Model.fitting_nc_amp_limit;
        if force_amp_scale_factor < force_scale_factor
            force_scale_factor = force_amp_scale_factor;
            limit_type(iMode) = "amp";
        end
    end

    initial_force_ratio(:,[2*iMode-1,2*iMode]) = force_ratio*force_scale_factor;

end

nc_mode_index = modes>1000;

if num_uncalibrated_modes > 0
    [r,~,~,E,sep_id,f] = Model.add_sep(initial_force_ratio,"lambda");
    
    removal_indicies = [];
    beyond_limit_index = E > Model.fitting_energy_limit;
    if ~isempty(Model.fitting_nc_amp_limit)
        beyond_amp_index = abs(f(nc_mode_index,:)) > Model.fitting_nc_amp_limit;
        beyond_limit_index = beyond_limit_index | beyond_amp_index;
    end
    num_seps = max(sep_id);
    for iSep = 1:num_seps
        sep_index = sep_id == iSep;
        beyond_sep = sep_index & beyond_limit_index;
        if nnz(beyond_sep) > 1
            removal_index = find(beyond_sep);
            removal_indicies = [removal_indicies,removal_index(2:end)]; %#ok<AGROW>
        end
    end
    r_all = r;
    f_all = f;
    E_all = E;
    sep_id_all = sep_id;

    r(:,removal_indicies) = [];
    f(:,removal_indicies)= [];
    E(:,removal_indicies) = [];
    sep_id(:,removal_indicies) = [];
    

    
end


for iMode = 1:num_uncalibrated_modes
    mode = uncalibrated_modes(iMode);
    mode_index = mode == modes;

    sep_span = (sep_id == 2*iMode-1 | sep_id == 2*iMode);


    eval_mode = Model.reduced_eigenvalues(mode_index);
  

    num_seps = 2;
    r_limit = zeros(1,num_seps);
    f_limit = zeros(1,num_seps);
    force_poly = cell(1,2);
    potential_poly = cell(1,2);
    for iSep = 1:num_seps
        sep_span = sep_id == (iSep+2*(iMode-1));
        E_sep = E(sep_span);
        r_sep = r(mode_index,sep_span);
        f_sep = f(mode_index,sep_span);

        switch limit_type(iMode)
            case "amp"
                f_diff = abs(f_sep) - abs(Model.nc_amplitude_limit);
                f_upper = f_diff;
                f_lower = f_diff;

                f_upper(f_upper < 0) = inf;
                f_lower(f_lower > 0) = -inf;

                [~,min_index] = max(f_lower);
                [~,max_index] = min(f_upper);
                bound_index = [min_index,max_index];
            case "energy"
                E_diff = E_sep - Model.energy_limit;
                E_upper = E_diff;
                E_lower = E_diff;

                E_upper(E_upper < 0) = inf;
                E_lower(E_lower > 0) = -inf;

                [~,min_index] = max(E_lower);
                [~,max_index] = min(E_upper);
                bound_index = [min_index,max_index];
        end

        num_loadcases = size(f_sep,2);
        max_force_degree = num_loadcases+1;
        force_degree = min(max_force_degree,11);
        Force_Poly_i = Polynomial(r_sep,f_sep,force_degree,"constraint",{"linear_force",eval_mode},"coupling","force","shift",1,"scale",1);

        if ~Model.Static_Options.follower_force
            Potential_Poly_i = integrate_polynomial(Force_Poly_i);
        else
            Potential_Poly_i = Polynomial(r_sep,E_sep,force_degree+1,"constraint",{"linear",[0,0]},"shift",1,"scale",1);
             % Potential_Poly_i = Polynomial(r_sep,E_sep,force_degree+1,"shift",1,"scale",1);
        end

        r_bound = r_sep(bound_index);
        r_interp = linspace(r_bound(1),r_bound(2));
        v_interp = Potential_Poly_i.evaluate_polynomial(r_interp);
        

        f_interp = Force_Poly_i.evaluate_polynomial(r_interp);

        switch limit_type(iMode)
            case "amp"
                [~,min_index] = min(abs(abs(f_interp)-abs(Model.nc_amplitude_limit)));
            case "energy"
                [~,min_index] = min(abs(v_interp - Model.energy_limit));
        end


        r_limit(1,iSep) = r_interp(min_index);
        f_limit(1,iSep) = f_interp(min_index);

        force_poly{iSep} = Force_Poly_i;
        potential_poly{iSep} = Potential_Poly_i;
    end


    Force_Calibration.force_limit{1,calibration_id}(iMode+num_calibrated_modes,:) = f_limit;
    Force_Calibration.calibrated_modes{1,calibration_id}(iMode+num_calibrated_modes,:) = mode;
    %---------------------------------------------------------%
    model_calibration_plot(mode,sep_id_all,iMode,r_all(mode_index,:),f_all(mode_index,:),f_limit,E_all,force_poly,potential_poly,Model)
    %---------------------------------------------------------%
end
Force_Calibration.Parameters = Model.Parameters;

save(GEOMETRY_PATH + "force_calibration","Force_Calibration")

calibrated_modes = Force_Calibration.calibrated_modes{1,calibration_id};
Model.calibrated_forces = zeros(num_modes,2);
% obj.calibrated_degree_limits = Force_Calibration.min_degree_data{1,calibration_id};

for iMode = 1:num_modes
    mode = modes(iMode);
    Model.calibrated_forces(iMode,:) = Force_Calibration.force_limit{1,calibration_id}(mode == calibrated_modes,:)*Calibration_Opts.force_overcalibration;
    % obj.calibrated_degree_limits{iMode}.force_applied_force = obj.calibrated_degree_limits{iMode}.force_applied_force./obj.calibrated_forces(iMode,:)';
    % obj.calibrated_degree_limits{iMode}.disp_applied_force = obj.calibrated_degree_limits{iMode}.disp_applied_force./obj.calibrated_forces(iMode,:)';
end
end