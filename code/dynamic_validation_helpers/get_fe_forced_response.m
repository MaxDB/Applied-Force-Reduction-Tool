function fe_forced_orbits = get_fe_forced_response(orbits,Rom,Force_Data,Damping_Data,Add_Output)
MIN_INC_SCALE_FACTOR = 1;
NUM_PERIODS = 50;

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
initial_time = zeros(1,num_orbits);
for iOrbit = 1:num_orbits
    orbit = orbits{iOrbit,1};
    t = orbit.tbp;
    r = orbit.xbp';
    period(iOrbit) = orbit.T;
    
    num_points = size(r,2);
    min_incs(iOrbit) = num_points*MIN_INC_SCALE_FACTOR;

    r_disp = r(1:num_modes,:);
    r_dot = r(num_modes + (1:num_modes),:);
    V = Rom.Potential_Polynomial.evaluate_polynomial(r_disp);
    [~,max_V_index] = max(V);

    % max_V_index = 1;
    initial_displacement = r_disp(:,max_V_index);
    initial_velocity = r_dot(:,max_V_index);
    initial_time(:,iOrbit) = t(max_V_index);

    initial_force(:,iOrbit) = Rom.Force_Polynomial.evaluate_polynomial(initial_displacement);

    physical_displacement(:,iOrbit) = Rom.expand(initial_displacement);
    physical_velocity(:,iOrbit) = Rom.expand_velocity(initial_displacement,initial_velocity);
end

switch Force_Data.type
    case "modal"
        FE_Force_Data.amplitude = Force_Data.amplitude;
        FE_Force_Data.harmonic_coefficients = [0,0,1];
        FE_Force_Data.alpha = Damping_Data.mass_factor;
        FE_Force_Data.beta = Damping_Data.stiffness_factor;
        FE_Force_Data.force_shape = Model.reduced_eigenvectors(:,Force_Data.mode_number);
end

[t_fom,x_fom,x_dot_fom,energy_fom] = Model.dynamic_simulation(physical_displacement,physical_velocity,initial_force,period,NUM_PERIODS,min_incs,initial_time,FE_Force_Data);

if ~iscell(x_fom)
    x_fom = {x_fom};
    t_fom = {t_fom};
    x_dot_fom = {x_dot_fom};
    energy_fom = {energy_fom};
end


frequency = zeros(1,num_orbits);
energy = zeros(1,num_orbits);
amplitude = zeros(num_modes,num_orbits);
periodicity = zeros(1,num_orbits);
switch Add_Output.type
    case "physical displacement"
        additional_dynamic_output = zeros(1,num_orbits);
end

r_transform = Model.reduced_eigenvectors'*Model.mass;
for iOrbit = 1:num_orbits
    t_fom_i = t_fom{iOrbit};
    x_fom_i = x_fom{iOrbit};
    x_dot_fom_i = x_dot_fom{iOrbit};
    energy_sim = energy_fom{iOrbit};
    
    orbit = orbits{iOrbit,1};
    period = orbit.T;
    time_range = [NUM_PERIODS-1,NUM_PERIODS]*period;
    in_range = t_fom_i >= time_range(1) & t_fom_i <= time_range(2);
    t_sim = t_fom_i(in_range);
    t_sim = t_sim - t_sim(1);
    x_sim = x_fom_i(:,in_range);
    x_dot_sim = x_dot_fom_i(:,in_range);

    potential_energy = energy_sim.potential(:,in_range(2:end));
    external_work = energy_sim.work(:,in_range(2:end));
    kinetic_energy = energy_sim.kinetic(:,in_range(2:end));
    dissipated_energy = energy_sim.dissipated(:,in_range(2:end));

    periodicity_error = norm(x_sim(:,end) - x_sim(:,1))/norm(x_sim(:,1));
    periodicity(1,iOrbit) = periodicity_error;
    % if ~converged
    %     warning("FE forced response " + iOrbit + " may not have converged")
    % end
        
    r_sim = r_transform*x_sim;

    frequency(1,iOrbit) = 2*pi/period;
    energy(1,iOrbit) = max(kinetic_energy+potential_energy);
    amplitude(:,iOrbit) = abs(max(r_sim,[],2) - min(r_sim,[],2))/2;

    switch Add_Output.type
        case "physical displacement"
            x_dof = x_sim(Add_Output.dof,:);
            additional_dynamic_output(:,iOrbit) = max(abs(x_dof),[],2);
    end
    % fe_forced_orbit.t = t_sim;
    % fe_forced_orbit.x = x_sim;
    % fe_forced_orbit.x_dot = x_dot_sim;
    % save()
end

fe_forced_orbits.frequency = frequency;
fe_forced_orbits.amplitude = amplitude;
fe_forced_orbits.energy = energy;
fe_forced_orbits.periodicity = periodicity;

switch Add_Output.type
    case "physical displacement"
        fe_forced_orbits.additional_dynamic_output = additional_dynamic_output;
end

% t_fom_1 = t_fom{1,1};
% x_fom_1 = x_fom{1,1};
% r_fom = Model.reduced_eigenvectors'*Model.mass*x_fom_1;
% 
% orbit_rom = orbits{1,1};
% t_rom = [];
% t_orbit = orbit_rom.tbp';
% orbit_period = t_orbit(end);
% t_orbit = t_orbit - initial_time;
% for iPeriod = 1:(NUM_PERIODS+1)
%     t_rom = [t_rom,t_orbit + orbit_period*(iPeriod-1)];
% end
% z_rom = repmat(orbit_rom.xbp',1,(NUM_PERIODS+1));
% r_rom = z_rom(1:num_modes,:);
% r_dot_rom = z_rom((1:num_modes) + 1,:);
% 
% figure
% hold on
% plot(t_rom,r_rom,"-k")
% plot(t_fom_1,r_fom,"--r")
% hold off


end



