function Static_Data = local_verification(Dyn_Data,sol_nums)
DELTA = 0.01;

%-------
Static_Data = load_static_data(Dyn_Data);

%--
collate_points_time_start = tic;
orbit_points = get_orbit_points(Dyn_Data,sol_nums);
total_points = size(orbit_points,2);
collate_points_time = toc(collate_points_time_start);
log_message = sprintf("%u trajectory points collated in %.1f seconds" ,total_points,collate_points_time);
logger(log_message,3)

select_points_time_start = tic;
orbit_points = select_orbit_points(orbit_points,DELTA);
verification_points = size(orbit_points,2);
select_points_time = toc(select_points_time_start);
log_message = sprintf("%u trajectory points selected in %.1f seconds" ,verification_points,select_points_time);
logger(log_message,3)




end


%----
function r_points = get_orbit_points(Dyn_Data,sol_nums)
num_sols = length(sol_nums);
num_modes = Dyn_Data.get_reduced_dimension;

r_points = zeros(num_modes,0);
for iSol = 1:num_sols
    sol_num = sol_nums(iSol);
    Sol = Dyn_Data.load_solution(sol_num);
    num_orbits = Sol.num_orbits;
    r_sol_points = zeros(num_modes,num_orbits*300);
    point_counter = 0;
    for iOrbit = 1:num_orbits
        Orbit = Dyn_Data.get_orbit(sol_num,iOrbit);
        r = Orbit.xbp(:,1:num_modes)';
        r_length = size(r,2);
        new_point_counter = point_counter + r_length;
        r_sol_points(:,(point_counter+1):new_point_counter) = r;
        point_counter = new_point_counter;
    end
    r_sol_points(:,(point_counter+1):end) = [];
    r_points = [r_points,r_sol_points]; %#ok<AGROW>
end
end
%----
function r_points = select_orbit_points(all_r_points,delta)

r_max = max(abs(all_r_points),[],2);
all_r_points = all_r_points./r_max;
%----
point_distance = vecnorm(all_r_points);
[sorted_point_distance,sort_index] = sort(point_distance,"ascend");
num_total_points = length(sorted_point_distance);

num_selected_points = 1;
r_points = zeros(size(all_r_points));
r_points_distance = zeros(size(sorted_point_distance));

for point_counter = 1:num_total_points
    selected_r_points = r_points(:,1:num_selected_points);
    selected_r_points_distance = r_points_distance(1:num_selected_points);
    
    point_index = point_counter;
    
    potential_r_point = all_r_points(:,sort_index(point_index));
    potential_r_point_distance = sorted_point_distance(point_index);
    
    comparison_index = abs(selected_r_points_distance - potential_r_point_distance) < delta;
    comp_distance = vecnorm(selected_r_points(:,comparison_index) - potential_r_point);
    if any(comp_distance < delta)
        continue
    end
    num_selected_points = num_selected_points + 1;
    r_points(:,num_selected_points) = potential_r_point;
    r_points_distance(num_selected_points) = potential_r_point_distance;
end
r_points(:,(num_selected_points+1):end) = [];
r_points = r_points.*r_max;
end