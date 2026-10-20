function Trajectory = get_full_order_trajectory(Model,duration,Nonconservative_Input,Simulation_Opts,Initial_Conditions)


period = 2*pi/Nonconservative_Input.frequency;
num_periods = duration/period;
min_incs = period/Simulation_Opts.max_time_step;
%---
if isempty(Initial_Conditions)
    t0 = 0;
    x0 = zeros(Model.num_dof,1);
    x_dot0 = zeros(Model.num_dof,1);
    f0 = zeros(Model.num_dof,1);
else
    t0 = Initial_Conditions.time;
    x0 = Initial_Conditions.disp;
    x_dot0 = Initial_Conditions.vel;
    f0 = Initial_Conditions.f0;
end

job_id = 1;


[t,x,x_dot,energy]  = Model.dynamic_simulation(x0,x_dot0,f0,period,num_periods,min_incs,t0,Nonconservative_Input,job_id);
Trajectory.t = t;
Trajectory.x = x;
Trajectory.x_dot = x_dot;
Trajectory.energy = energy;
end

