function num_jobs = get_current_parallel_jobs
    parallel_pool = gcp("nocreate");
    if isempty(parallel_pool)
        num_jobs = 0;
    else
        num_jobs = parallel_pool.NumWorkers;
    end
end