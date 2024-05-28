function job_force_ratios = split_abaqus_jobs(force_ratios,num_loadcases,Static_Opts)
max_parallel_jobs = Static_Opts.max_parallel_jobs;
mininum_job_loadcases = Static_Opts.minimum_job_loadcases;

num_seps = size(force_ratios,2);
maximum_jobs = ceil(num_seps*num_loadcases/mininum_job_loadcases);
num_parallel_jobs = min([maximum_jobs,max_parallel_jobs,num_seps]);

job_force_ratios = cell(1,num_parallel_jobs);

remaining_seps = num_seps;
remaining_groups = num_parallel_jobs;
last_sep = 0;
for iJob = 1:num_parallel_jobs
    group_size = ceil(remaining_seps/remaining_groups);
    next_seps = last_sep + (1:group_size);
    job_force_ratios{1,iJob} = force_ratios(:,next_seps);
    
    last_sep = next_seps(end);
    remaining_seps = remaining_seps - group_size;
    remaining_groups = remaining_groups - 1;
end
end