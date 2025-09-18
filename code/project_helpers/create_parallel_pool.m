function parallel_pool = create_parallel_pool(num_workers)
resources = "Processes";
parpool_input = {resources,num_workers};

max_workers = 512;
if isstring(num_workers)
    switch num_workers
        case "current"
            return
        case "max"
            num_workers = max_workers;
            parpool_input(2) = [];
        otherwise
            error("Unknown option: '"+ num_workers + "'")
    end
end

parpool_input(end+1) = {"SpmdEnabled"};
parpool_input(end+1) = {false};

parallel_pool = gcp("nocreate");

%create pool if it doesn't exist
if isempty(parallel_pool) && num_workers > 1
    parallel_pool = parpool(parpool_input{:});
    return
end

%don't create pool if it doesn't exist and there are not multiple workers
if isempty(parallel_pool) && num_workers <= 1 
    return
end

%should be unreachable
if isempty(parallel_pool)
    error("")
end

if num_workers == parallel_pool.NumWorkers
    return
end

%delete pool if it does exist and there are not multiple workers
if num_workers <= 1
    delete(parallel_pool)
    parallel_pool = gcp("nocreate");
    return
end


if num_workers == max_workers
    profile_workers = parallel_pool.Cluster.NumWorkers;
    if profile_workers == parallel_pool.NumWorkers
        return
    end
end

if num_workers ~= parallel_pool.NumWorkers
    delete(parallel_pool)
    parallel_pool = parpool(parpool_input{:});
    return
end

error("Unimplemented case")
end


