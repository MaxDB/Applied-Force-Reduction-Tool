function get_full_order_forced_response(Model,Nonconservative_Input,solution_num,reference_sol)
NUM_SIM_PERIODS = 200;
%----
num_dof = size(Nonconservative_Input.force_shape,1);
x_0 = zeros(num_dof,1);
x_dot_0 = zeros(num_dof,1);
f_0 = zeros(num_dof,1);


min_incs = 100;

%----
frequency = Nonconservative_Input.frequency;
periods = 2*pi./frequency;
%----

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
for iOrbit = 1:num_orbits
    orbit_start_time = tic;
    period = periods(iOrbit);
    
    if ~isempty(reference_sol)
        if isfield(reference_sol,"orbit_subset")
            orbit_id = reference_sol.orbit_subset(iOrbit);
        else
            orbit_id = iOrbit;
        end
        ref_orbit = Dyn_Data.get_orbit(reference_sol.sol_num,orbit_id);
        num_modes = size(ref_orbit.xbp,2)/2;
        r_orbit = ref_orbit.xbp(:,1:num_modes)';
        r_dot_orbit = ref_orbit.xbp(:,(1:num_modes) + num_modes)';
        x_0 = Dyn_Data.Dynamic_Model.expand(r_orbit(:,1));
        x_dot_0 = Dyn_Data.Dynamic_Model.expand_velocity(r_orbit(:,1),r_dot_orbit(:,1));
        fr_0 = Dyn_Data.Dynamic_Model.Force_Polynomial.evaluate_polynomial(r_orbit(:,1));
        f_0 = Model.mass*Model.reduced_eigenvectors*fr_0;
    end

    initial_time = 0;
    job_id = [1,1];
    num_periods = NUM_SIM_PERIODS;
    
    Nonconservative_Input.fe_output = "none";


    [t,x,~,energy]  = Model.dynamic_simulation(x_0,x_dot_0,f_0,period,num_periods,min_incs,initial_time,Nonconservative_Input,job_id);
    
    job_id = [1,2];
    Nonconservative_Input.fe_output = "all";
    num_periods = 1.5;
    [t_per,x_per,x_dot_per,energy_per]  = Model.dynamic_simulation([],[],[],period,num_periods,min_incs,t(end),Nonconservative_Input,job_id);
    
    %-- test
    % figure; plot(t(2:end),energy.potential+energy.kinetic)
    % hold on
    % plot([t(end),t_per],[energy.potential(end)+energy.kinetic(end),energy_per.potential+energy_per.kinetic])
    %--

    % check for convergence
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

    save(data_path+"\sol"+iOrbit,"Orbit_Data")
    orbit_time = toc(orbit_start_time);
    log_message = sprintf("FOM FR %.3f: %u/%u in %.1f seconds" ,periodicity_error,iOrbit,num_orbits,orbit_time);
    logger(log_message,2)
end

end