function get_full_order_forced_response(Model,Nonconservative_Input,solution_num,reference_sol,varargin)
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

num_sim_periods = [];
for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "num_periods"
            num_sim_periods = keyword_values{arg_counter};
        case "num_workers"
            num_parallel_workers = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------%
current_pool = gcp("nocreate");
if ~isempty(current_pool) && current_pool.NumWorkers ~= num_parallel_workers
    delete(current_pool)
    parpool(num_parallel_workers);
elseif isempty(current_pool)
    parpool(num_parallel_workers);
end

%----
min_incs = 100;

%----
frequency = Nonconservative_Input.frequency;
periods = 2*pi./frequency;
%----
if isempty(num_sim_periods)
    quality_factor = sqrt(Model.reduced_eigenvalues(1))/Nonconservative_Input.alpha;
    num_sim_periods = ceil(100*quality_factor);
end
%--
data_path = split(Model.get_data_path,"\");
data_path = join(data_path(1:2),"\") + "\dynamic_sol_" + solution_num;
if ~isfolder(data_path)
    mkdir(data_path)
end
%---
if ~isempty(reference_sol)
    Dyn_Data = initalise_dynamic_data(reference_sol.name);
    Model = Dyn_Data.Dynamic_Model.Model;
end

num_orbits = length(periods);
num_parallel_groups = ceil(num_orbits/num_parallel_workers);
num_dof = Model.num_dof;

for iGroup = 1:num_parallel_groups
    group_start_time = tic;
    orbit_group = ((iGroup-1)*num_parallel_workers+1):(iGroup*num_parallel_workers);
    orbit_group(orbit_group>num_orbits) = [];
    num_group_orbits = length(orbit_group);
    x_0 = zeros(num_dof,num_group_orbits);
    x_dot_0 = zeros(num_dof,num_group_orbits);
    f_0 = zeros(num_dof,num_group_orbits);

    group_ids = zeros(1,num_group_orbits);
    for iOrbit = 1:num_group_orbits
        orbit_id = orbit_group(iOrbit);
        period = periods(orbit_id);

        if ~isempty(reference_sol)
            if isfield(reference_sol,"orbit_subset")
                orbit_id = reference_sol.orbit_subset(orbit_id);
            end
            ref_orbit = Dyn_Data.get_orbit(reference_sol.sol_num,orbit_id);
            num_modes = size(ref_orbit.xbp,2)/2;
            r_orbit = ref_orbit.xbp(:,1:num_modes)';
            r_dot_orbit = ref_orbit.xbp(:,(1:num_modes) + num_modes)';
            x_0(:,iOrbit) = Dyn_Data.Dynamic_Model.expand(r_orbit(:,1));
            x_dot_0(:,iOrbit) = Dyn_Data.Dynamic_Model.expand_velocity(r_orbit(:,1),r_dot_orbit(:,1));
            fr_0 = Dyn_Data.Dynamic_Model.Force_Polynomial.evaluate_polynomial(r_orbit(:,1));
            f_0(:,iOrbit) = Model.mass*Model.reduced_eigenvectors*fr_0;
        end
        group_ids(iOrbit) = orbit_id;
    end
    initial_time = 0;
    job_id = [(1:num_group_orbits)',ones(num_group_orbits,1)];
    num_periods = num_sim_periods;

    Nonconservative_Input.fe_output = "none";
    sim_periods = period*ones(1,num_group_orbits);
    sim_num_periods = num_periods*ones(1,num_group_orbits);
    sim_min_incs = min_incs*ones(1,num_group_orbits);
    sim_initial_time = initial_time*ones(1,num_group_orbits);
    [t,x,~,energy]  = Model.dynamic_simulation(x_0,x_dot_0,f_0,sim_periods,sim_num_periods,sim_min_incs,sim_initial_time,Nonconservative_Input,job_id);

    job_id(:,2) = job_id(:,2) + 1;
    Nonconservative_Input.fe_output = "all";
    if iscell(t)
        sim_t_start = cellfun(@(t_cell) t_cell(end),t);
    else
        sim_t_start = t(end);
    end
    num_periods = 1.5;
    sim_num_periods = num_periods*ones(1,num_group_orbits);
    [t_per_sim,x_per_sim,x_dot_per_sim,energy_per_sim]  = Model.dynamic_simulation(zeros(0,num_group_orbits),zeros(0,num_group_orbits),zeros(0,num_group_orbits),sim_periods,sim_num_periods,sim_min_incs,sim_t_start,Nonconservative_Input,job_id);

    %-- test
    % figure; plot(t(2:end),energy.potential+energy.kinetic)
    % hold on
    % plot([t(end),t_per],[energy.potential(end)+energy.kinetic(end),energy_per.potential+energy_per.kinetic])
    %--
    for iOrbit = 1:num_group_orbits
        if iscell(t_per_sim)
            t_per = t_per_sim{iOrbit};
            x_per = x_per_sim{iOrbit};
            x_dot_per = x_dot_per_sim{iOrbit};
            energy_per = energy_per_sim{iOrbit};
        else
            t_per = t_per_sim;
            x_per = x_dot_per_sim;
            x_dot_per = x_dot_per_sim;
            energy_per_sim = energy_per;
        end
        % check for convergence
        if isempty(t_per)
            log_message = sprintf("Abaqus error: orbit id %u",group_ids(iOrbit));
            logger(log_message,2)
            continue
        end
        final_time = t_per(end);
        period_start = final_time-period;

        t_index = find(t_per >= period_start);
        if t_index(1) > 1
            t_index = [t_index(1)-1,t_index];
        end
        t_periodic = t_per(t_index);
        x_periodic = x_per(:,t_index);
        x_dot_periodic = x_dot_per(:,t_index);
        potential_periodic = energy_per.potential(:,t_index);
        kinetic_periodic = energy_per.kinetic(:,t_index);


        for iDof = 1:num_dof
            x_periodic(iDof,1) = interp1(t_periodic(1:2),x_periodic(iDof,1:2),period_start);
            x_dot_periodic(iDof,1) = interp1(t_periodic(1:2),x_dot_periodic(iDof,1:2),period_start);
        end
        potential_periodic(1) = interp1(t_periodic(1:2),potential_periodic(1:2),period_start);
        kinetic_periodic(1) = interp1(t_periodic(1:2),kinetic_periodic(1:2),period_start);

        t_periodic(1) = period_start;

        x_start = x_periodic(:,1);
        x_end = x_periodic(:,end);



        periodicity_error = norm(x_end-x_start)/norm(x_start);

        converged = periodicity_error < 1e-3;

        if ~converged
            warning("FRC not converged")
        end

        %KE test
        % num_time_points = size(t_periodic,2);
        % kinetic_energy = zeros(1,num_time_points);
        % for iTime = 1:num_time_points
        %     kinetic_energy(iTime) = 0.5*x_dot_periodic(:,iTime)'*Model.mass*x_dot_periodic(:,iTime);
        % end

        %

        Orbit_Data.time = t_periodic;
        Orbit_Data.disp = x_periodic;
        Orbit_Data.vel = x_dot_periodic;
        Orbit_Data.energy = potential_periodic + kinetic_periodic;
        Orbit_Data.periodicity_error = periodicity_error;
        Orbit_Data.period = period;

        save(data_path+"\sol"+orbit_group(iOrbit),"Orbit_Data")
        orbit_time = toc(group_start_time);
        log_message = sprintf("FOM FR %.3f: %u/%u in %.1f seconds" ,periodicity_error,orbit_group(iOrbit),num_orbits,orbit_time);
        logger(log_message,2)
    end
end

end