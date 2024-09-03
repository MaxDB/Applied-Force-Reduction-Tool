function static_dataset_verificiation_plot(Static_Data)
PLOT_LEVEL = 1;
NUM_DISP_PLOTS = 4;

load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end

Rom = Reduced_System(Static_Data);

ax = plot_static_data("energy",Static_Data);
Rom.Potential_Polynomial.plot_polynomial(ax);

ax = plot_static_data("force",Static_Data);
Rom.Force_Polynomial.plot_polynomial(ax);

num_dofs = Static_Data.Model.num_dof;
disp_outputs = randi(num_dofs,[NUM_DISP_PLOTS,1]);

ax = plot_static_data("displacement",Static_Data,"outputs",disp_outputs);
Rom.Physical_Displacement_Polynomial.plot_polynomial(ax,disp_outputs);

end