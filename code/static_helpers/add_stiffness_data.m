function Static_Data = add_stiffness_data(Static_Data,Static_Opts)
Static_Opts.additional_data = "stiffness";

Static_Data.Model = Static_Data.Model.update_static_opts(Static_Opts);
Static_Data.additional_data_type = "stiffness";

x_disp = Static_Data.physical_displacement;
force = Static_Data.restoring_force;

Point_Data.initial_disp = x_disp;
Point_Data.initial_force = force;
[~,~,~,~,stiffness] = Static_Data.Model.add_point([],"stiffness",Point_Data);

Static_Data = Static_Data.update_additional_data(stiffness);
end