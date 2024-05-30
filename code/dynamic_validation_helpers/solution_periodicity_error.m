function periodicity_error = solution_periodicity_error(orbit,Rom)
MIN_INC_SCALE_FACTOR = 1;
Model = Rom.Model;
period = orbit.T;
t = orbit.tbp;
r = orbit.xbp';



num_modes = size(r,1)/2;
num_points = size(r,2);
min_incs = num_points*MIN_INC_SCALE_FACTOR;

r_dot = r(num_modes + (1:num_modes),:);
approx_ke = zeros(1,num_points);
for iPoint = 1:num_points
    approx_ke(iPoint) = r_dot(:,iPoint)'*r_dot(:,iPoint);
end
[~,min_ke_index] = min(approx_ke);


initial_displacement = r(1:num_modes,min_ke_index);
initial_velocity = r(num_modes + (1:num_modes),min_ke_index);

initial_force = Rom.Force_Polynomial.evaluate_polynomial(initial_displacement);

physical_displacement = Rom.expand(initial_displacement);
physical_velocity = Rom.expand_velocity(initial_displacement,initial_velocity);


[t_fom,x_fom,x_dot_fom] = Model.dynamic_simulation(physical_displacement,physical_velocity,initial_force,period,min_incs);
periodicity_error = norm(x_fom(:,end) - x_fom(:,1))/norm(x_fom(:,1));
end