function fe_forced_orbits = get_fe_forced_response(orbits,Rom,Force_Data,Damping_Data,Add_Output)
MIN_INC_SCALE_FACTOR = 1;
NUM_PERIODS = 100;
MAX_PERIODICITY_ERROR = 1e-10; %convergence criteria
MAX_ITERATIONS = 1000;
MIN_ITERATIONS = 2; %give it some time for stuff to go wrong

DEBUG_PLOT = 1;

Model = Rom.Model;
% num_dofs = Model.num_dof;
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
        FE_Force_Data.force_shape = Model.mass*Model.reduced_eigenvectors(:,Force_Data.mode_number);
    case "point force"
        FE_Force_Data.amplitude = Force_Data.amplitude;
        FE_Force_Data.harmonic_coefficients = [0,0,1];
        FE_Force_Data.alpha = Damping_Data.mass_factor;
        FE_Force_Data.beta = Damping_Data.stiffness_factor;

        num_dofs = Model.num_dof;
        dof_map = zeros(num_dofs,1);
        dof_map(Model.node_mapping(:,1) == Force_Data.dof) = 1;
        FE_Force_Data.force_shape = dof_map;
end

max_parallel_jobs = Model.Static_Options.max_parallel_jobs;
num_parallel_jobs = min(max_parallel_jobs,num_orbits);
orbit_groups = split_abaqus_jobs(1:num_orbits,1,max_parallel_jobs,1);
Model.Static_Options.max_parallel_jobs = 1;

r_transform = Model.reduced_eigenvectors'*Model.mass;

frequency_group = cell(1,num_parallel_jobs);
energy_group = cell(1,num_parallel_jobs);
amplitude_group = cell(1,num_parallel_jobs);
periodicity_group = cell(1,num_parallel_jobs);
simulated_periods_group = cell(1,num_parallel_jobs);
switch Add_Output.type
    case "physical displacement"
        additional_dynamic_output_group = cell(1,num_parallel_jobs);
        node_map = Model.node_mapping;
        monitored_dof = node_map(node_map(:,1) == Add_Output.dof,2);
end

reset_temp_directory()

% parfor (iJob = 1:num_parallel_jobs,max_parallel_jobs)
for iJob = 1:num_parallel_jobs
    orbit_group = orbit_groups{iJob};
    num_group_orbits = size(orbit_group,2);

    frequency = zeros(1,num_group_orbits);
    energy = zeros(1,num_group_orbits);
    amplitude = zeros(num_modes,num_group_orbits);
    periodicity = zeros(1,num_group_orbits);
    additional_dynamic_output = zeros(1,num_group_orbits);
    simulated_periods = zeros(1,num_group_orbits);

    for iOrbit = 1:num_group_orbits
        orbit_sim_time_start = tic;

        orbit_id = orbit_group(iOrbit);
        orbit = orbits{orbit_id,1};
        t = orbit.tbp;
        r = orbit.xbp';
        period = orbit.T;

        num_points = size(r,2);
        min_incs = num_points*MIN_INC_SCALE_FACTOR;

        r_disp = r(1:num_modes,:);
        r_dot = r(num_modes + (1:num_modes),:);
        % V = Rom.Potential_Polynomial.evaluate_polynomial(r_disp);
        % [~,start_index] = max(V);
        start_index = 1;

        % max_V_index = 1;
        initial_displacement = r_disp(:,start_index);
        initial_velocity = r_dot(:,start_index);
        initial_time = t(start_index);

        initial_force = Rom.Force_Polynomial.evaluate_polynomial(initial_displacement);

        physical_displacement = Rom.expand(initial_displacement);
        physical_velocity = Rom.expand_velocity(initial_displacement,initial_velocity);
        
        if DEBUG_PLOT
            fig_all = figure; %#ok<*UNRCH>
            fig_all.Name = "Job " + iJob + ", orbit " + iOrbit;
            tiledlayout("flow")
            ax_all = cell(1,num_modes);
            for iMode = 1:num_modes
                ax_all{iMode} = nexttile;
                box(ax_all{iMode},"on")
                xlabel(ax_all{iMode},"t (s)")
                ylabel(ax_all{iMode},"q_{"+Model.reduced_modes(iMode) + "}")
                hold(ax_all{iMode},"on")
            end

            switch Add_Output.type
                case "physical displacement"
                    fig_physical = figure;
                    fig_physical.Name = "Job " + iJob + ", orbit " + iOrbit;
                    tiledlayout("flow")
                    ax_physical = nexttile;
                    hold(ax_physical,"on")
                    box(ax_physical,"on")
                    xlabel(ax_physical,"t (s)")
                    ylabel(ax_physical,"x_{"+Add_Output.dof+"}")
            end
        end
        
        x0 = physical_displacement;
        for iStep = 1:MAX_ITERATIONS
            orbit_sim_step_start = tic;

            job_id = [iJob,iStep];
            [t_fom,x_fom,~,energy_fom] = Model.dynamic_simulation(physical_displacement,physical_velocity,initial_force,period,NUM_PERIODS,min_incs,initial_time,FE_Force_Data,job_id);
            
           
            delete("temp\dynamic_analysis_" + iJob + "_" + (iStep - 2)+".*")

            x_sim = x_fom;
            xN = x_fom(:,end);
            periodicity_error = norm(xN - x0)/norm(x0);
            x0 = xN;

            if DEBUG_PLOT
                r_fom = r_transform*x_fom;
                for iMode = 1:num_modes
                    plot(ax_all{iMode},t_fom+period*NUM_PERIODS*(iStep-1),r_fom(iMode,:))
                end
                
                switch Add_Output.type
                    case "physical displacement"
                        x_dof = x_fom(monitored_dof,:);
                        plot(ax_physical,t_fom+period*NUM_PERIODS*(iStep-1),x_dof)
                end
                % t_sim_norm = t_fom(:,period_span(1):period_span(2))/period;
                % t_sim_norm = t_sim_norm - t_sim_norm(1);
                % r_sim = r_transform*x_sim;
                % [t_shift,r_shift] = shift_orbit(t_sim_norm,r_sim);
                % for iMode = 1:num_modes
                %     plot(ax_period{iMode},t_shift,r_shift(iMode,:))
                % end

                drawnow
                file_name = "temp\J" + iJob + "O" + iOrbit + "_";
                saveas(fig_all,file_name+"all.fig")
                saveas(fig_physical,file_name+"period.fig")
            end

            converged = periodicity_error < MAX_PERIODICITY_ERROR;

            orbit_sim_step_time = toc(orbit_sim_step_start);
            log_message = sprintf("Job " + iJob + ", orbit " + iOrbit + ": %i/%i periods in %.1f seconds with %.3f periodicity error" ,NUM_PERIODS,iStep*NUM_PERIODS,orbit_sim_step_time,periodicity_error);
            logger(log_message,2)

            if converged && iStep >= MIN_ITERATIONS
                break
            end
        end

        potential_energy = energy_fom.potential;
        kinetic_energy = energy_fom.kinetic;

        % external_work = energy_sim.work(:,in_range(2:end));
        % dissipated_energy = energy_sim.dissipated(:,in_range(2:end));

        r_sim = r_transform*x_sim;
        simulated_periods(1,iOrbit) = iStep*NUM_PERIODS;
        periodicity(1,iOrbit) = periodicity_error;
        frequency(1,iOrbit) = 2*pi/period;
        energy(1,iOrbit) = max(kinetic_energy+potential_energy);
        amplitude(:,iOrbit) = abs(max(r_sim,[],2) - min(r_sim,[],2))/2;

        switch Add_Output.type
            case "physical displacement"
                x_dof = x_sim(monitored_dof,:);
                additional_dynamic_output(:,iOrbit) = max(abs(x_dof),[],2);
        end

        orbit_sim_time = toc(orbit_sim_time_start);
        log_message = sprintf("Orbit " + iOrbit + ": %i period simulation complete in %.1f seconds" ,NUM_PERIODS*iStep,orbit_sim_time);
        logger(log_message,1)
    end

    frequency_group{1,iJob} = frequency;
    amplitude_group{1,iJob} = amplitude;
    energy_group{1,iJob} = energy;
    periodicity_group{1,iJob} = periodicity;
    additional_dynamic_output_group{1,iJob} = additional_dynamic_output;
    simulated_periods_group{1,iJob} = simulated_periods;
end


fe_forced_orbits.frequency = [frequency_group{:}];
fe_forced_orbits.amplitude = [amplitude_group{:}];
fe_forced_orbits.energy = [energy_group{:}];
fe_forced_orbits.periodicity = [periodicity_group{:}];
fe_forced_orbits.simulated_periods = [simulated_periods_group{:}];

switch Add_Output.type
    case "physical displacement"
        fe_forced_orbits.additional_dynamic_output = [additional_dynamic_output_group{:}];
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



