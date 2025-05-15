function max_parallel_jobs = maximum_viable_jobs(num_seps,num_sep_input_lines,additional_job_input_lines,Static_Opts)
max_inp_lines = Static_Opts.max_lines_per_job;
if num_sep_input_lines > max_inp_lines
    warning("System likely too big for effective parallelisation (assuming 32GB of RAM). If working with more memory increase the max_lines_per_job setting in static_opts")
end

max_parallel_jobs = Static_Opts.max_parallel_jobs;
for iJob = 1:max_parallel_jobs
    concurrent_lines = iJob*ceil(num_seps/iJob)*num_sep_input_lines + iJob*additional_job_input_lines;
    if concurrent_lines > max_inp_lines
        max_parallel_jobs = iJob-1;
        break
    end
end


end