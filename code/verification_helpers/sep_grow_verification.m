function Static_Data = sep_grow_verification(Static_Data)
% add set of points at degree boundary
% add points until ROM has converged up to that bondary
% move to next boundary or actually energy limit
start_degree_index = 2;
sep_density = 5;

Model = Static_Data.Model;
min_degree_data = Model.calibrated_degree_limits;

num_r_modes = size(Model.reduced_modes,2);
max_degree = 3;

for iMode = 1:num_r_modes
    Mode_Degree_Data = min_degree_data{iMode};
    force_degree = Mode_Degree_Data.force_degrees;
    disp_degree = Mode_Degree_Data.disp_degrees;
    max_mode_degree = max(max(force_degree),max(disp_degree));
    max_degree = max(max_degree,max_mode_degree);
end

num_steps = size(force_degree,2) + 1;
step_scale_factor = zeros(1,num_steps-1) + inf;
step_energy_limit = zeros(1,num_steps-1) + inf;

for iMode = 1:num_r_modes
    Mode_Degree_Data = min_degree_data{iMode};
    force_af = mean(Mode_Degree_Data.force_applied_force,1);
    disp_af = mean(Mode_Degree_Data.disp_applied_force,1);

    force_energy = Mode_Degree_Data.force_degree_energy;
    disp_energy = Mode_Degree_Data.disp_degree_energy;

    energy_limit = min(force_energy,disp_energy);
    lim_index = (disp_energy < force_energy);
    applied_force = [force_af;disp_af];
    applied_force = applied_force([~lim_index;lim_index])';
    
    step_scale_factor = min(step_scale_factor,applied_force);
    step_energy_limit = min(step_energy_limit,energy_limit);
end
step_scale_factor = [step_scale_factor,1];
step_energy_limit = [step_energy_limit,Model.energy_limit];
step_fitting_energy_limit = step_energy_limit*Model.Calibration_Options.energy_overfit;


for iStep = 1:num_steps
    % set model
    Sub_Model = Model;
    Sub_Model.energy_limit = step_energy_limit(iStep);
    Sub_Model.fitting_energy_limit = step_fitting_energy_limit(iStep);
    Static_Data.Model = Sub_Model;
    

    % create scaffold
    Static_Data = Static_Data.create_scaffold("density",sep_density,"sf",step_scale_factor(iStep));
    % plot_static_data("energy",Static_Data,"plot seps",0)

    % verify limited stess manifold
    
    Static_Data = sep_from_origin_verification(Static_Data);
end
end