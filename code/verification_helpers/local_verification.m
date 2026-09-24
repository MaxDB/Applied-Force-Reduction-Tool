function Static_Data = local_verification(Dyn_Data,sol_nums,varargin)
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

Verification_Opts = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "verification_opts"
            Verification_Opts = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------
DELTA = 0.01;

%-------
Static_Data = load_static_data(Dyn_Data);
if ~isempty(Verification_Opts)
    Static_Data = Static_Data.update_verification_opts(Verification_Opts);
end
%--
%get points
collate_points_time_start = tic;
orbit_points = get_orbit_points(Dyn_Data,sol_nums);
total_points = size(orbit_points,2);
collate_points_time = toc(collate_points_time_start);

log_message = sprintf("%u trajectory points collated in %.1f seconds" ,total_points,collate_points_time);
logger(log_message,3)

%select point subset
select_points_time_start = tic;
[orbit_points,r_max] = select_orbit_points(orbit_points,DELTA);
num_verification_points = size(orbit_points,2);
select_points_time = toc(select_points_time_start);

log_message = sprintf("%u trajectory points selected in %.1f seconds" ,num_verification_points,select_points_time);
logger(log_message,3)

%verify points
error_time_start = tic;
[interpolation_error,Rom_One] = get_interpolation_error(orbit_points,Static_Data);
max_error = max(interpolation_error);
num_failures = nnz(interpolation_error > 1);
error_time = toc(error_time_start);

log_message = sprintf("%u failed points with maximum error %.1f found in %.1f seconds" ,num_failures,max_error,error_time);
logger(log_message,3)


%select loadcases to add
point_add_time_start = tic;
added_points = select_points_to_add(orbit_points,interpolation_error,r_max,DELTA);
new_loads = Rom_One.Force_Polynomial.evaluate_polynomial(added_points);
[r_new,x_new,force_new,energy_new,additional_data_new] = Rom_One.Model.add_point(new_loads,Static_Data.additional_data_type);
Static_Data = Static_Data.update_data(r_new,x_new,force_new,energy_new,0,additional_data_new);
point_add_time = toc(point_add_time_start);

log_message = sprintf("%u points added in %.1f seconds" ,size(added_points,2),point_add_time);
logger(log_message,3)

%recheck error
error_time_start = tic;
interpolation_error = get_interpolation_error(orbit_points,Static_Data);
max_error = max(interpolation_error);
num_failures = nnz(interpolation_error > 1);
error_time = toc(error_time_start);

log_message = sprintf("%u failed points with maximum error %.1f found in %.1f seconds" ,num_failures,max_error,error_time);
logger(log_message,3)

%clean up
system_name = get_system_name(Static_Data);
delete_dynamic_data(system_name);
Static_Data.save_data;

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
function [r_points,r_max] = select_orbit_points(all_r_points,delta)

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
r_points(:,1) = []; %remove origin
r_points = r_points.*r_max;
end

%-----------------

function [interpolation_error,Rom_One] = get_interpolation_error(orbit_points,Static_Data)
current_degree = Static_Data.verified_degree;
Verification_Opts = Static_Data.Verification_Options;

Rom_One = Reduced_System(Static_Data,"degree",current_degree);
Rom_Two = Reduced_System(Static_Data,"degree",current_degree+2);
[force_error,disp_error] = verify_points(orbit_points,Rom_One,Rom_Two);
norm_force_error = force_error/Verification_Opts.maximum_interpolation_error(1);
norm_disp_error = disp_error/Verification_Opts.maximum_interpolation_error(2);
interpolation_error = max(norm_force_error,norm_disp_error);

end

%-----------------

function added_points = select_points_to_add(orbit_points,interpolation_error,r_max,delta)
error_indicies = find(interpolation_error>1);
if isempty(error_indicies)
    added_points = [];
    return
end
[~,indices_order] = sort(interpolation_error(error_indicies),"descend");

error_points = orbit_points(:,error_indicies(indices_order));
scaled_error_points = error_points./r_max;

num_points = size(scaled_error_points,2);
added_points = zeros(size(scaled_error_points));
added_points(:,1) = scaled_error_points(:,1);
num_added_points = 1;
for iPoint = 2:num_points
    scaled_error_point = scaled_error_points(:,iPoint); 
    point_distance = vecnorm(added_points(:,1:num_added_points) - scaled_error_point);
    if any(point_distance < 10*delta)
        continue
    end
    num_added_points = num_added_points + 1;
    added_points(:,num_added_points) = scaled_error_point;
end
added_points(:,(num_added_points+1):end) = [];
added_points = added_points.*r_max;
end