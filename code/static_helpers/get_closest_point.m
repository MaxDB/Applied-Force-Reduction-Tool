function Closest_Points = get_closest_point(Static_Data,new_loads)
%not very sophisticated but search is based on number of points and force dimension so should be alright
restoring_force = Static_Data.get_dataset_values("restoring_force");
displacement = Static_Data.get_dataset_values("physical_displacement");

num_new_loads = size(new_loads,2);
initial_disp = zeros(size(displacement,1),num_new_loads);
initial_force = zeros(size(restoring_force,1),num_new_loads);
for iLoad = 1:num_new_loads
    restoring_force_diff = restoring_force - new_loads(:,iLoad);
    force_distance = sum(restoring_force_diff.^2,1);
    [~,closest_point_id] = min(force_distance);
    initial_disp(:,iLoad) = displacement(:,closest_point_id);
    initial_force(:,iLoad) = restoring_force(:,closest_point_id);
end

Closest_Points.initial_disp = initial_disp;
Closest_Points.initial_force = initial_force;
end