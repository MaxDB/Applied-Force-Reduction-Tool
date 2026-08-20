function static_dataset_verification_plot(Static_Data)
PLOT_LEVEL = 1;
NUM_DISP_PLOTS = 4;

load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end

Rom = Reduced_System(Static_Data);
if length(Rom.Model.reduced_modes) > 2
    return
end
limit_data = {Rom.Potential_Polynomial,Rom.Model.energy_limit,Rom.Model.nc_amplitude_limit};
num_dofs = Static_Data.Model.num_dof;

ax = plot_static_data("energy",Static_Data);
Rom.Potential_Polynomial.plot_polynomial("axes",ax,"potential",limit_data);


if size(Rom.Force_Polynomial,1) > 10
    force_outputs = randi(num_dofs,[NUM_DISP_PLOTS,1]);
    ax = plot_static_data("force",Static_Data,"outputs",force_outputs);
    Rom.Force_Polynomial.plot_polynomial("axes",ax,"potential",limit_data,"outputs",force_outputs);
else
    ax = plot_static_data("force",Static_Data);
    Rom.Force_Polynomial.plot_polynomial("axes",ax,"potential",limit_data);
end

disp_outputs = randi(num_dofs,[NUM_DISP_PLOTS,1]);

ax = plot_static_data("displacement",Static_Data,"outputs",disp_outputs);
Rom.Physical_Displacement_Polynomial.plot_polynomial("axes",ax,"potential",limit_data,"outputs",disp_outputs);

end