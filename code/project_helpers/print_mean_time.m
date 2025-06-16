function print_mean_time(time,name,exclude_data)
if nargin == 2
    exclude_data = [];
end
time(exclude_data) = [];

num_time_points = size(time,2);
num_outputs = size(time,1);

mean_time = mean(time,2);
max_time_diff = max(time,[],2) - mean_time;
min_time_diff = min(time,[],2) - mean_time;

if num_outputs == 1
fprintf("%-24s: %.1f +%.1f -%.1f (s) over %i points \n",name,mean_time,max_time_diff,abs(min_time_diff),num_time_points);
else
    for iOutput = 1:num_outputs
        fprintf("%-20s (%i): %.1f +%.1f -%.1f (s) over %i points \n",name,iOutput,mean_time(iOutput),max_time_diff(iOutput),abs(min_time_diff(iOutput)),num_time_points);
    end
end