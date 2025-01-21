classdef Validated_Backbone_Solution 
    properties
        validation_modes
        Validation_Options

        h_amplitude
        corrected_low_modal_amplitude
        low_modal_amplitude
        h_energy
        h_modal_energy_fraction
        h_stability
        h_force_amplitude
        r_force_amplitude
        h_abs_mean
        r_abs_mean
    end

    methods
        function obj = Validated_Backbone_Solution(Rom,BB_Sol,Validated_BB_Settings)
            
            Validation_Opts.validation_algorithm = "h_frequency";
            Validation_Opts.get_stability = 1;
            obj = obj.update_validation_opts(Validation_Opts);
            
            solution_num = Validated_BB_Settings.solution_num;
            L_modes = Validated_BB_Settings.L_modes;

            validation_start = tic;
            validation_rom_start = tic;
            

            load_static_data_start = tic;
            Static_Data = load_static_data(Rom);

            load_static_data_time = toc(load_static_data_start);
            log_message = sprintf("Static dataset loaded: %.1f seconds" ,load_static_data_time);
            logger(log_message,3)

            Static_Data = Static_Data.add_validation_data(L_modes);
            % save(data_path + "Static_Data.mat","Static_Data","-v7.3")

            L_modes = Static_Data.Dynamic_Validation_Data.current_L_modes;

            obj.validation_modes = L_modes;
            obj = obj.preallocate_analysis_outputs(BB_Sol);
            Validation_Rom = Reduced_System(Static_Data);

            validation_rom_time = toc(validation_rom_start);
            log_message = sprintf("Validation ROM created: %.1f seconds" ,validation_rom_time);
            logger(log_message,2)

            validation_solution_start = tic;
            obj = solve_h_prediction(obj,BB_Sol,Validation_Rom,solution_num);
            % switch H_ALGORITHM
            %     case "harmonic_balance"
            %         obj = h_harmonic_balance(BB_Sol,Validation_Rom,solution_num);
            %     case "time_solution"
            %         obj = h_time_solution(obj,BB_Sol,Validation_Rom,solution_num);
            % end
            validation_solution_time = toc(validation_solution_start);
            log_message = sprintf("Validation equations solved: %.1f seconds" ,validation_solution_time);
            logger(log_message,2)

            validation_time = toc(validation_start);
            log_message = sprintf("Solution validated: %.1f seconds" ,validation_time);
            logger(log_message,1)
        end
        %-----------------------------------------------------------------%
        function obj = preallocate_analysis_outputs(obj,BB_Sol)
            num_orbits = BB_Sol.num_orbits;
            
            num_validation_modes = size(obj.validation_modes,2);
            num_r_modes = size(BB_Sol.amplitude,1);

            num_h_modes = num_validation_modes + num_r_modes;

            obj.h_amplitude = zeros(num_h_modes,num_orbits);
            obj.corrected_low_modal_amplitude = zeros(num_h_modes,num_orbits);
            obj.low_modal_amplitude = zeros(num_h_modes,num_orbits);
            obj.h_energy = zeros(1,num_orbits);
            obj.h_stability = zeros(1,num_orbits);
            obj.h_force_amplitude = zeros(num_h_modes,num_orbits);
            obj.r_force_amplitude = zeros(num_h_modes,num_orbits);
            obj.h_abs_mean = zeros(num_h_modes,num_orbits);
            obj.r_abs_mean = zeros(num_h_modes,num_orbits);

        end
        %-----------------------------------------------------------------%
        function obj = analyse_h_solution(obj,Displacment,Velocity,Eom_Terms,orbit_stability,Validation_Analysis_Inputs,orbit_num)
            
            get_amplitude = @(x) (max(x,[],2) - min(x,[],2))/2;
            %--
            r = Displacment.r;
            h = Displacment.h;

            r_dot = Velocity.r_dot;
            h_dot = Velocity.h_dot;
            
            % h_stiff = Eom_Terms.h_stiffness;
            h_force = Eom_Terms.h_force;
            %--
            h_amp = get_amplitude(h);
            
            %--
            R_Disp_Poly = Validation_Analysis_Inputs.Physical_Disp_Data.Poly;
            x_r = R_Disp_Poly.evaluate_polynomial(r);

            L_disp_transform = Validation_Analysis_Inputs.L_disp_transform;
            g_l = L_disp_transform*x_r;

            g_h = [r;g_l];
            g_amp = get_amplitude(g_h);

            %--
            q_l = h+g_h;
            q_amp = get_amplitude(q_l);

            %--
            % h_ke_approx = 0.5*1
            % [~,energy_point_span] = min(abs(h_dot))
            energy_point_span = 1:size(h,2);
            num_points = size(energy_point_span,2);

            r_energy = r(:,energy_point_span);
            r_dot_energy = r_dot(:,energy_point_span);
            h_disp_energy = h(:,energy_point_span);
            h_dot_energy = h_dot(:,energy_point_span);

            num_r_modes = size(r_energy,1);
            num_h_modes = size(h,1);
            
             
            r_force = Validation_Analysis_Inputs.Force_Poly.evaluate_polynomial(r_energy);
            
            h_potential = zeros(1,num_points);
            for iPoint = 1:num_points
                h_i = h(:,iPoint);
                r_force_hat_i = [r_force(:,iPoint);zeros(num_h_modes-num_r_modes,1)];
                r_force_gradient_i = Validation_Analysis_Inputs.H_Force_Poly.evaluate_polynomial(r_energy(:,iPoint));
                h_potential(iPoint) = (r_force_hat_i + r_force_gradient_i*h_i)'*h_i;
            end


            potential_tilde = Validation_Analysis_Inputs.Potential_Poly.evaluate_polynomial(r_energy);
            potential_hat = h_potential + potential_tilde;
            [max_potential,energy_index] = max(potential_hat); 
            %--
            
            r_force_amp = get_amplitude(r_force);
            r_force_amp((num_r_modes+1):num_h_modes) = 0;
            % validation_force(1:num_r_modes,:) = validation_force(1:num_r_modes,:) + r_force;
            % h_force_amp = get_amplitude(validation_force);
            %--
            [~,max_ke_hat] = h_kinetic_energy(r_energy(:,energy_index),r_dot_energy(:,energy_index),h_disp_energy(:,energy_index),h_dot_energy(:,energy_index),Validation_Analysis_Inputs);
            % %--
            % energy_hat = ke_hat + potential_hat;
            %energy_hat = mean(energy_hat)
            energy_hat = max_potential + max_ke_hat;
                
            %--
            % TEST
            % h_energy = get_h_energy(r_energy,r_dot_energy,h_disp_energy,h_dot_energy,Validation_Analysis_Inputs);
            %---
            % h_force(h_force<1e-6) = 0;
            h_zero_force = get_amplitude(h_force);
            h_force_amp = h_zero_force;
            % r_force_amp((num_r_modes+1):num_h_modes) = get_amplitude(arrayfun(@norm,r_force));
            %--
            mean_abs = @(x) mean(abs(x),2);
            mean_abs_h = mean_abs(h);
            mean_abs_r2 = mean_abs(g_h);

            %--
            obj.h_amplitude(:,orbit_num) = h_amp;
            obj.corrected_low_modal_amplitude(:,orbit_num) = g_amp;
            obj.low_modal_amplitude(:,orbit_num) = q_amp;
            obj.h_energy(:,orbit_num) = energy_hat;
            obj.h_stability(:,orbit_num) = orbit_stability;
            obj.h_force_amplitude(:,orbit_num) = h_force_amp;
            obj.r_force_amplitude(:,orbit_num) = r_force_amp;
            obj.h_abs_mean(:,orbit_num) = mean_abs_h;
            obj.r_abs_mean(:,orbit_num) = mean_abs_r2;
        end
        %-----------------------------------------------------------------%
        function obj = update_validation_opts(obj,Validation_Opts)
            Default_Continuation_Opts = read_default_options("validation");
            obj.Validation_Options = update_options(Default_Continuation_Opts,obj.Validation_Options,Validation_Opts);
        end
    end

end