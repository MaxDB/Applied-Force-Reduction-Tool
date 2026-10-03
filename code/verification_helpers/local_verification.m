function [Static_Data,Static_Time_Data] = local_verification(Dyn_Data,sol_nums,varargin)
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

Verification_Opts = [];
Static_Time_Data.static_time = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "verification_opts"
            Verification_Opts = keyword_values{arg_counter};
        case "static_time"
            Static_Time_Data = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------
DELTA = 0.01;
MAXIMUM_DEGREE = 13;
MAX_POINTS = 0.25;
max_loadcases = 200;
% figure;
% ax = axes;
% hold(ax,"on")
% box(ax,"on")
% xlabel(ax,"r_1")
% ylabel(ax,"r_2")
%-------
Static_Data = load_static_data(Dyn_Data);
if ~isempty(Verification_Opts)
    Static_Data = Static_Data.update_verification_opts(Verification_Opts);
end
%--
if isempty(Static_Time_Data.static_time)
    max_points = ceil(size(Static_Data,2)*MAX_POINTS);
else
    points_fit = Static_Time_Data.max_loadcases;
    time_fit = Static_Time_Data.static_time;
    if isscalar(time_fit)
        points_fit = [0,points_fit];
        time_fit = [0,time_fit];
    end
    coeffs = lsqminnorm([time_fit;ones(size(points_fit))]',points_fit');
    max_points = 3*ceil(Static_Time_Data.dynamic_time(end)*coeffs(1) + coeffs(2));

    max_points = max(max_points,ceil(size(Static_Data,2)*MAX_POINTS));
end 
max_points = min(max_points,max_loadcases);
%--
select_points_time_start = tic;
%get points
collate_points_time_start = tic;
orbit_points = get_orbit_points(Dyn_Data,sol_nums);
total_points = size(orbit_points,2);
collate_points_time = toc(collate_points_time_start);

log_message = sprintf("%u trajectory points collated in %.1f seconds" ,total_points,collate_points_time);
logger(log_message,3)

%select point subset
[orbit_points,r_max] = select_orbit_points(orbit_points,DELTA);
num_verification_points = size(orbit_points,2);
select_points_time = toc(select_points_time_start);

log_message = sprintf("%u trajectory points selected in %.1f seconds" ,num_verification_points,select_points_time);
logger(log_message,2)

% plot(ax,orbit_points(1,:),orbit_points(2,:),"kx")

error_time_start = tic;
num_r_modes = Static_Data.get_reduced_dimension;
num_dataset_points = size(Static_Data,2);
max_force_degree = get_max_poly_degree("force",num_r_modes,num_dataset_points,MAXIMUM_DEGREE);
max_disp_degree = get_max_poly_degree("displacement",num_r_modes,num_dataset_points,MAXIMUM_DEGREE);
% 
% force_degree = max(Static_Data.verified_degree(1)-2,3);
% disp_degree = max(Static_Data.verified_degree(2)-2,3);

force_degree = Static_Data.verified_degree(1);
disp_degree = Static_Data.verified_degree(2);

num_degree_pairs = min(max_force_degree+2 - force_degree,max_disp_degree+2 - disp_degree)/2;
maximum_force_pair_errors = zeros(1,num_degree_pairs);
maximum_disp_pair_errors = zeros(1,num_degree_pairs);
interpolation_errors = zeros(num_verification_points,num_degree_pairs);
degree_pairs = zeros(2,num_degree_pairs);

for iPair = 1:num_degree_pairs
    %verify points
   
    error_iteration_time_start = tic;
    log_message = sprintf("Comparing %s and %s degree force polynomials" , ...
        ordinal_suffix(force_degree),ordinal_suffix(force_degree+2));
    logger(log_message,4)
    log_message = sprintf("Comparing %s and %s degree displacement polynomials" , ...
        ordinal_suffix(disp_degree),ordinal_suffix(disp_degree+2));
    logger(log_message,4)



    Static_Data.verified_degree = [force_degree,disp_degree];
    [interpolation_error,max_force_error,max_disp_error] = get_interpolation_error(orbit_points,Static_Data);
    error_time = toc(error_time_start);
    
 

    maximum_force_pair_errors(iPair) = max_force_error;
    maximum_disp_pair_errors(iPair) = max_disp_error;
    
    interpolation_errors(:,iPair) = interpolation_error;
    degree_pairs(:,iPair) = [force_degree;disp_degree];
    
    if max_force_error > 1
    force_degree = force_degree + 2;
    end

    if max_disp_error > 1
    disp_degree = disp_degree + 2;
    end

    error_iteration_time = toc(error_iteration_time_start);
    log_message = sprintf("Max force error: %.2f and max disp error: %.2f in %.1f seconds" ,max_force_error,max_disp_error,error_iteration_time);
    logger(log_message,4)

    if max_force_error < 1 && max_disp_error < 1
        remove_index = (iPair + 1):num_degree_pairs;
        maximum_force_pair_errors(remove_index) = [];
        maximum_disp_pair_errors(remove_index) = [];
        interpolation_errors(:,remove_index) = [];
        degree_pairs(:,remove_index) = [];
        break
    end
end
% plot(ax,orbit_points(1,interpolation_error>1),orbit_points(2,interpolation_error>1),"rx")

max_error = max(interpolation_errors,[],1);
[~,min_error_pair] = min(max_error);
force_degree = degree_pairs(1,min_error_pair);
disp_degree = degree_pairs(1,min_error_pair);
interpolation_error = interpolation_errors(:,min_error_pair);
num_failures = nnz(interpolation_error > 1);

log_message = sprintf("%u failed points with maximum error %.1f found in %.1f seconds" ,num_failures,max_error(min_error_pair),error_time);
logger(log_message,3)

if num_failures == 0
    Static_Data.verified_degree = [force_degree,disp_degree];
    Static_Data.save_data;
    return
end

%select loadcases to add
point_add_time_start = tic;
Rom_One = Reduced_System(Static_Data,"degree",[force_degree,disp_degree]);
existing_points = Static_Data.reduced_displacement;
added_points = select_points_to_add(orbit_points,interpolation_error,r_max,DELTA,max_points,existing_points);
new_loads = Rom_One.Force_Polynomial.evaluate_polynomial(added_points);


local_verification_plot(Static_Data,orbit_points,interpolation_error,added_points)


% plot(ax,added_points(1,:),added_points(2,:),"bx")

%add loadcases
static_time_start = tic;
[r_new,x_new,force_new,energy_new,additional_data_new] = Rom_One.Model.add_point(new_loads,Static_Data.additional_data_type,[]);
Static_Data = Static_Data.update_data(r_new,x_new,force_new,energy_new,0,additional_data_new);
static_time = toc(static_time_start);
point_add_time = toc(point_add_time_start);

num_added_points = size(added_points,2);
log_message = sprintf("%u points added in %.1f seconds" ,num_added_points,point_add_time);
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

num_iteration = length(Static_Time_Data.dynamic_time);
Static_Time_Data.static_time(num_iteration) = static_time;
Static_Time_Data.max_loadcases(num_iteration) = num_added_points;
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

function [interpolation_error,max_force_error,max_disp_error] = get_interpolation_error(orbit_points,Static_Data)
current_degree = Static_Data.verified_degree;
Verification_Opts = Static_Data.Verification_Options;

Rom_One = Reduced_System(Static_Data,"degree",current_degree);
Rom_Two = Reduced_System(Static_Data,"degree",current_degree+2);
[force_error,disp_error] = verify_points(orbit_points,Rom_One,Rom_Two);
norm_force_error = force_error/Verification_Opts.maximum_interpolation_error(1);
norm_disp_error = disp_error/Verification_Opts.maximum_interpolation_error(2);
interpolation_error = max(norm_force_error,norm_disp_error);
max_force_error = max(norm_force_error);
max_disp_error = max(norm_disp_error);
end

%-----------------

function added_points = select_points_to_add(orbit_points,interpolation_error,r_max,delta,max_points,existing_points)
distance_reduction = 0.5;
min_distance = 10*delta;

error_indicies = find(interpolation_error>1);
if isempty(error_indicies)
    added_points = [];
    return
end
[~,indices_order] = sort(interpolation_error(error_indicies),"descend");

error_points = orbit_points(:,error_indicies(indices_order));
num_error_points = size(error_points,2);
if num_error_points < max_points
    added_points = error_points;
    return
end

scaled_error_points = error_points./r_max;



num_existing_points = size(existing_points,2);
scaled_existing_points = existing_points./r_max; 

added_points = zeros(size(scaled_error_points));
added_points = [scaled_existing_points,added_points];
num_added_points = num_existing_points + 1;
max_points = max_points + num_existing_points;

added_points(:,num_added_points) = scaled_error_points(:,1);
max_points_added = false;
inc_counter = 0;
while ~max_points_added
    num_points = size(scaled_error_points,2);
    removal_index = false(1,size(scaled_error_points,2));
    for iPoint = 1:num_points
        scaled_error_point = scaled_error_points(:,iPoint);
        point_distance = vecnorm(added_points(:,1:num_added_points) - scaled_error_point);
        if any(point_distance < min_distance)
            continue
        end

        removal_index(iPoint) = true;
        num_added_points = num_added_points + 1;
        added_points(:,num_added_points) = scaled_error_point;

        if num_added_points == max_points
            max_points_added = true;
            break
        end
    end
    if num_added_points - max_points < 16
        min_distance = min_distance*distance_reduction;
    else
        max_points_added = true;
    end
    scaled_error_points(:,removal_index) = [];
    inc_counter = inc_counter + 1;
    if inc_counter == 3
        max_points_added = true;
    end
end

added_points(:,(num_added_points+1):end) = [];
added_points(:,1:num_existing_points) = [];
added_points = added_points.*r_max;
end