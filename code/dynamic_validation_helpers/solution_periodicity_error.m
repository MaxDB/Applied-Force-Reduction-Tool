function periodicity_error = solution_periodicity_error(orbits,Rom)
MIN_INC_SCALE_FACTOR = 1;
NUM_PERIODS = 1;

Model = Rom.Model;
num_dofs = Model.num_dof;
num_modes = length(Model.reduced_modes);

num_orbits = size(orbits,1);
if num_orbits == 1
    orbits = {orbits};
end

physical_displacement = zeros(num_dofs,num_orbits);
physical_velocity = zeros(num_dofs,num_orbits);
initial_force = zeros(num_modes,num_orbits);
period = zeros(1,num_orbits);
min_incs = zeros(1,num_orbits);

for iOrbit = 1:num_orbits
    orbit = orbits{iOrbit,1};
    period(iOrbit) = orbit.T;
    t = orbit.tbp;
    r = orbit.xbp';
    
    num_points = size(r,2);
    min_incs(iOrbit) = num_points*MIN_INC_SCALE_FACTOR;

    r_dot = r(num_modes + (1:num_modes),:);
    approx_ke = zeros(1,num_points);
    for iPoint = 1:num_points
        approx_ke(iPoint) = r_dot(:,iPoint)'*r_dot(:,iPoint);
    end
    [~,min_ke_index] = min(approx_ke);


    initial_displacement = r(1:num_modes,min_ke_index);
    initial_velocity = r(num_modes + (1:num_modes),min_ke_index);

    initial_force(:,iOrbit) = Rom.Force_Polynomial.evaluate_polynomial(initial_displacement);

    physical_displacement(:,iOrbit) = Rom.expand(initial_displacement);
    physical_velocity(:,iOrbit) = Rom.expand_velocity(initial_displacement,initial_velocity);
end

[t_fom,x_fom,x_dot_fom] = Model.dynamic_simulation(physical_displacement,physical_velocity,initial_force,period,NUM_PERIODS,min_incs);

if ~iscell(x_fom)
    x_fom = {x_fom};
end

periodicity_error = zeros(1,num_orbits);
for iOrbit = 1:num_orbits
    disp_fom = x_fom{1,iOrbit};
    periodicity_error(1,iOrbit) = norm(disp_fom(:,end) - disp_fom(:,1))/norm(disp_fom(:,1));

end

end