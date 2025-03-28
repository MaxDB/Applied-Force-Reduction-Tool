function validation_dataset_verificiation_plot(Static_Data)
PLOT_LEVEL = 2;
NUM_OUTPUT_PLOTS = 4;


load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end

num_dofs = Static_Data.Model.num_dof;

h_modes = Static_Data.get_current_h_data;
num_h_modes = size(h_modes,2);

Rom = Reduced_System(Static_Data);

stiffness_outputs = randi(num_h_modes,[NUM_OUTPUT_PLOTS,2]);
ax = plot_static_data("h_stiffness",Static_Data,"outputs",stiffness_outputs);
Rom.Low_Frequency_Stiffness_Polynomial.plot_polynomial("axes",ax,"outputs",stiffness_outputs);

disp_output_rows = randi(num_dofs,[NUM_OUTPUT_PLOTS,1]);
disp_output_cols = randi(num_h_modes,[NUM_OUTPUT_PLOTS,1]);
disp_outputs = [disp_output_rows,disp_output_cols];
ax = plot_static_data("h_displacement_gradient",Static_Data,"outputs",disp_outputs);
Rom.Low_Frequency_Coupling_Gradient_Polynomial.plot_polynomial("axes",ax,"outputs",disp_outputs);

% ax = plot_static_data("energy",Static_Data);
% Rom.Potential_Polynomial.plot_polynomial(ax);
% 
% ax = plot_static_data("force",Static_Data);
% Rom.Force_Polynomial.plot_polynomial(ax);
% 
% 
% disp_outputs = randi(num_dofs,[1,NUM_DISP_PLOTS]);
% 
% ax = plot_static_data("displacement",Static_Data,"outputs",disp_outputs);
% Rom.Physical_Displacement_Polynomial.plot_polynomial(ax,disp_outputs);

end