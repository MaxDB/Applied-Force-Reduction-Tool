classdef Validated_Backbone_Solution 
    properties
        validation_modes

        h_amplitude
        corrected_low_modal_amplitude
        low_modal_amplitude
        h_energy
        h_modal_energy_fraction
        h_stability
        h_force_amplitude
        r_force_amplitude
    end

    methods
        function obj = Validated_Backbone_Solution(Rom,BB_Sol,Validated_BB_Settings)
            H_ALGORITHM = "time_solution";
            
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
            switch H_ALGORITHM
                case "harmonic_balance"
                    obj = h_harmonic_balance(BB_Sol,Validation_Rom,solution_num);
                case "time_solution"
                    obj = h_time_solution(obj,BB_Sol,Validation_Rom,solution_num);
            end
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

        end
        %-----------------------------------------------------------------%
        function obj = analyse_h_solution(obj,r,r_dot,h,h_dot,orbit_stability,Validation_Analysis_Inputs,orbit_num)
            
            get_amplitude = @(x) (max(x,[],2) - min(x,[],2))/2;
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
            energy_point_span = 1:size(h,2);
            num_points = size(energy_point_span,2);

            r_energy = r(:,energy_point_span);
            r_dot_energy = r_dot(:,energy_point_span);
            h_disp_energy = h(:,energy_point_span);
            h_dot_energy = h_dot(:,energy_point_span);

            H_Stiff_Poly = Validation_Analysis_Inputs.H_Force_Poly;
            num_r_modes = size(r_energy,1);
            num_h_modes = size(h,1);
            
            h_potential = zeros(1,num_points);
            h_force = zeros(num_h_modes,num_points);
            for iPoint = 1:num_points
                h_stiffness = H_Stiff_Poly.evaluate_polynomial(r_energy(:,iPoint));
                h_force(:,iPoint) = h_stiffness*h_disp_energy(:,iPoint);
                h_potential(:,iPoint) = 0.5*h_force(:,iPoint)'*h_disp_energy(:,iPoint);
            end

            potential_tilde = Validation_Analysis_Inputs.Potential_Poly.evaluate_polynomial(r_energy);
            potential_hat = h_potential + potential_tilde;
            %--
            r_force = Validation_Analysis_Inputs.Force_Poly.evaluate_polynomial(r_energy);
            r_force_amp = get_amplitude(r_force);
            r_force_amp((num_r_modes+1):num_h_modes) = 0;
            h_force(1:num_r_modes,:) = h_force(1:num_r_modes,:) + r_force;
            h_force_amp = get_amplitude(h_force);
            %--
            [ke_tilde,ke_hat] = h_kinetic_energy(r_energy,r_dot_energy,h_disp_energy,h_dot_energy,Validation_Analysis_Inputs);
            %--
            energy_hat = ke_hat + potential_hat;
                

            %--
            obj.h_amplitude(:,orbit_num) = h_amp;
            obj.corrected_low_modal_amplitude(:,orbit_num) = g_amp;
            obj.low_modal_amplitude(:,orbit_num) = q_amp;
            obj.h_energy(:,orbit_num) = mean(energy_hat);
            obj.h_stability(:,orbit_num) = orbit_stability;
            obj.h_force_amplitude(:,orbit_num) = h_force_amp;
            obj.r_force_amplitude(:,orbit_num) = r_force_amp;
        end
        %-----------------------------------------------------------------%
    end

end