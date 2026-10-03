function local_verification_plot(Static_Data,orbit_points,interpolation_error,added_points)
PLOT_LEVEL = 2;

load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end

reduced_dims = size(orbit_points,1);
if reduced_dims > 3
    return
end

error_points = orbit_points(:,interpolation_error>1);
orbit_points(:,interpolation_error>1) = [];
error_points(:,ismember(error_points',added_points',"rows")) = [];



point_style = {"Marker",".","Color","k","LineStyle","none","MarkerSize",4};
error_style = {"Marker",".","Color",get_plot_colours(3),"LineStyle","none","MarkerSize",6};
added_point_style = {"Marker","x","Color",get_plot_colours(4),"LineStyle","none","MarkerSize",6,"LineWidth",2};

figure
tiledlayout
ax = nexttile;
box(ax,"on")
hold(ax,"on")
switch reduced_dims
    case 2
        plot(ax,orbit_points(1,:),orbit_points(2,:),point_style{:})
        plot(ax,error_points(1,:),error_points(2,:),error_style{:})
        plot(ax,added_points(1,:),added_points(2,:),added_point_style{:})
    case 3
        plot3(ax,orbit_points(1,:),orbit_points(2,:),orbit_points(3,:),point_style{:})
        plot3(ax,error_points(1,:),error_points(2,:),error_points(3,:),error_style{:})
        plot3(ax,added_points(1,:),added_points(2,:),added_points(3,:),added_point_style{:})
        zlabel(ax,"r_3")
end
hold(ax,"off")
xlabel(ax,"r_1")
ylabel(ax,"r_2")


%-------
ax = nexttile;
box(ax,"on")
hold(ax,"on")
reduced_disp = Static_Data.reduced_displacement;
switch reduced_dims
    case 2
        plot(ax,reduced_disp(1,:),reduced_disp(2,:),point_style{:})
        plot(ax,added_points(1,:),added_points(2,:),added_point_style{:})
    case 3
        plot3(ax,reduced_disp(1,:),reduced_disp(2,:),reduced_disp(3,:),point_style{:})
        plot3(ax,added_points(1,:),added_points(2,:),added_points(3,:),added_point_style{:})
        zlabel(ax,"r_3")
end
hold(ax,"off")
xlabel(ax,"r_1")
ylabel(ax,"r_2")
end