function get_full_order_forced_response(Model,Nonconservative_Input,solution_num)

%----
num_dof = size(Nonconservative_Input.force_shape,1);
x_0 = zeros(num_dof,1);
x_dot_0 = zeros(num_dof,1);
f_0 = zeros(num_dof,1);
initial_time = 0;
job_id = 1;

min_incs = 100;
num_periods = 100;
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


num_orbits = length(periods);
for iOrbit = 1:num_orbits
    period = periods(iOrbit);
    [t,x,x_dot,energy]  = Model.dynamic_simulation(x_0,x_dot_0,f_0,period,num_periods,min_incs,initial_time,Nonconservative_Input,job_id);

    % check for convergence
    final_time = t(end);
    period_start = final_time-period;

    t_index = find(t >= period_start);
    if t_index(1) > 1
        t_index = [t_index(1)-1,t_index];
    end
    t_periodic = t(t_index);
    x_periodic = x(:,t_index);

    x_start = zeros(num_dof,1);
    for iDof = 1:num_dof
        x_start(iDof,:) = interp1(t_periodic,x_periodic(iDof,:),period_start);
    end
    x_end = x_periodic(:,end);

    periodicity_error = norm(x_end-x_start)/norm(x_start);

    converged = periodicity_error < 1e-3;

    if ~converged
        warning("FRC not converged")
    end

    energy_index = t_index;
    if size(x,2) > size(energy.potential,2)
        energy_index = energy_index - 1;
    end

    Orbit_Data.time = t_periodic;
    Orbit_Data.disp = x_periodic;
    Orbit_Data.energy = energy.potential(energy_index) + energy.kinetic(energy_index);
    Orbit_Data.periodicity_error = periodicity_error;
    Orbit_Data.period = period;

    save(data_path+"\sol"+iOrbit,"Orbit_Data")
end

end