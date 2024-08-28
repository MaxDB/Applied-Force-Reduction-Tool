function fe_forced_orbits = get_fe_forced_response(orbits,Rom,Force_Data,Damping_Data,Add_Output)
MIN_INC_SCALE_FACTOR = 1;
NUM_PERIODS = 10;
MAX_PERIODICITY_ERROR = 1e-3;
MAX_ITERATIONS = 50;

Model = Rom.Model;
num_dofs = Model.num_dof;
num_modes = length(Model.reduced_modes);

num_orbits = size(orbits,1);
if num_orbits == 1
    orbits = {orbits};
end


switch Force_Data.type
    case "modal"
        FE_Force_Data.amplitude = Force_Data.amplitude;
        FE_Force_Data.harmonic_coefficients = [0,0,1];
        FE_Force_Data.alpha = Damping_Data.mass_factor;
        FE_Force_Data.beta = Damping_Data.stiffness_factor;
        FE_Force_Data.force_shape = Model.reduced_eigenvectors(:,Force_Data.mode_number);
end

max_parallel_jobs = Model.Static_Options.max_parallel_jobs;
orbit_groups = split_abaqus_jobs(1:num_orbits,1,max_parallel_jobs,1);
Model.Static_Options.max_parallel_jobs = 1;

r_transform = Model.reduced_eigenvectors'*Model.mass;

frequency = zeros(1,num_orbits);
energy = zeros(1,num_orbits);
amplitude = zeros(num_modes,num_orbits);
periodicity = zeros(1,num_orbits);
switch Add_Output.type
    case "physical displacement"
        additional_dynamic_output = zeros(1,num_orbits);
end
% parfor (iOrbit = 1:max_parallel_jobs,Static_Opts.max_parallel_jobs)
for iJob = 1:min(max_parallel_jobs,num_orbits)
    orbit_group = orbit_groups{iJob};
    num_group_orbits = size(orbit_group,2);
    for iOrbit = 1:num_group_orbits
        orbit_id = orbit_group(iOrbit);
        orbit = orbits{orbit_id,1};
        t = orbit.tbp;
        r = orbit.xbp';
        period = orbit.T;

        num_points = size(r,2);
        min_incs = num_points*MIN_INC_SCALE_FACTOR;

        r_disp = r(1:num_modes,:);
        r_dot = r(num_modes + (1:num_modes),:);
        V = Rom.Potential_Polynomial.evaluate_polynomial(r_disp);
        [~,max_V_index] = max(V);

        % max_V_index = 1;
        initial_displacement = r_disp(:,max_V_index);
        initial_velocity = r_dot(:,max_V_index);
        initial_time = t(max_V_index);

        initial_force = Rom.Force_Polynomial.evaluate_polynomial(initial_displacement);

        physical_displacement = Rom.expand(initial_displacement);
        physical_velocity = Rom.expand_velocity(initial_displacement,initial_velocity);
        
        figure
        tiledlayout("flow")
        for iMode = 1:num_modes
            ax{iMode} = nexttile;
            hold(ax{iMode},"on")
        end
        for iStep = 1:MAX_ITERATIONS

            [t_fom,x_fom,x_dot_fom,energy_fom] = Model.dynamic_simulation(physical_displacement,physical_velocity,initial_force,period,NUM_PERIODS,min_incs,initial_time,FE_Force_Data,iJob);
            
            time_range = [NUM_PERIODS-1,NUM_PERIODS]*period;
            in_range = t_fom >= time_range(1) & t_fom <= time_range(2);
            % t_sim = t_fom_i(in_range);
            % t_sim = t_sim - t_sim(1);
            x_sim = x_fom(:,in_range);
            % x_dot_sim = x_dot_fom_i(:,in_range);
            

            r_fom = r_transform*x_fom;
            for iMode = 1:num_modes
                plot(ax{iMode},t_fom+period*NUM_PERIODS*(iStep-1),r_fom(iMode,:))
            end

            periodicity_error = norm(x_sim(:,end) - x_sim(:,1))/norm(x_sim(:,1));

            converged = periodicity_error < MAX_PERIODICITY_ERROR;
            if converged
                break
            end

            physical_displacement = x_sim(:,end);
            physical_velocity = x_dot_fom(:,end);
            
            initial_displacement = r_transform*physical_displacement;
            
            initial_force = Rom.Force_Polynomial.evaluate_polynomial(initial_displacement);
        end

        potential_energy = energy_fom.potential(:,in_range(2:end));
        kinetic_energy = energy_fom.kinetic(:,in_range(2:end));

        % external_work = energy_sim.work(:,in_range(2:end));
        % dissipated_energy = energy_sim.dissipated(:,in_range(2:end));

        r_sim = r_transform*x_sim;
        
        periodicity(1,iOrbit) = periodicity_error;
        frequency(1,iOrbit) = 2*pi/period;
        energy(1,iOrbit) = max(kinetic_energy+potential_energy);
        amplitude(:,iOrbit) = abs(max(r_sim,[],2) - min(r_sim,[],2))/2;

        switch Add_Output.type
            case "physical displacement"
                x_dof = x_sim(Add_Output.dof,:);
                additional_dynamic_output(:,iOrbit) = max(abs(x_dof),[],2);
        end
    end
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



