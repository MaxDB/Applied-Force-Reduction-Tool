function orbit_groups = split_orbit_jobs(num_orbits,num_jobs)
    orbit_groups = zeros(num_jobs,2);
    remaining_jobs = num_jobs;
    first_orbit = 1;
    for iJob = 1:num_jobs
        orbits_per_job = ceil(num_orbits/remaining_jobs);
        orbit_groups(iJob,:) = [first_orbit,first_orbit + orbits_per_job - 1];
        num_orbits = num_orbits - orbits_per_job;
        remaining_jobs = remaining_jobs - 1;
        first_orbit = first_orbit + orbits_per_job;
    end
end