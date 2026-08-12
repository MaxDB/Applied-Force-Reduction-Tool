function validation_dataset_verificiation_plot(Static_Data)
PLOT_LEVEL = 2;
NUM_OUTPUT_PLOTS = 4;
view_angles = length(Static_Data.Model.reduced_modes) + 1;

plotting_level = load_log_level("plot");
if plotting_level < PLOT_LEVEL
    return
end

num_dofs = Static_Data.Model.num_dof;

h_modes = Static_Data.get_current_h_data;
num_h_modes = size(h_modes,2);

Rom = Reduced_System(Static_Data);

if num_h_modes^2 <= NUM_OUTPUT_PLOTS
    mode_index = 1:num_h_modes;
    stiffness_outputs = zeros(NUM_OUTPUT_PLOTS,2);
    stiffness_outputs(:,1) = repelem(mode_index,2);
    stiffness_outputs(:,2) = repmat(mode_index',2,1);
else
    stiffness_outputs = randi(num_h_modes,[NUM_OUTPUT_PLOTS,2]);
end
ax = plot_static_data("h_stiffness",Static_Data,"outputs",stiffness_outputs);
Rom.Low_Frequency_Stiffness_Polynomial.plot_polynomial("axes",ax,"outputs",stiffness_outputs);
for iAx = 1:length(ax)
    view(ax{iAx},view_angles);
end

disp_output_rows = randi(num_dofs,[NUM_OUTPUT_PLOTS,1]);
disp_output_cols = randi(num_h_modes,[NUM_OUTPUT_PLOTS,1]);
disp_outputs = [disp_output_rows,disp_output_cols];
ax = plot_static_data("h_displacement_gradient",Static_Data,"outputs",disp_outputs);
Rom.Low_Frequency_Coupling_Gradient_Polynomial.plot_polynomial("axes",ax,"outputs",disp_outputs);
for iAx = 1:length(ax)
    view(ax{iAx},view_angles);
end


end