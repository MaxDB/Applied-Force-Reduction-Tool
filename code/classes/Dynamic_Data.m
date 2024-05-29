classdef Dynamic_Data
    %Creates and stores dynamic results: periodic orbits, backbone curves, FRFs ect... 
    properties
        Dynamic_Model
        Continuation_Options

        num_solutions
        solution_types
        
        solution_labels
        frequency
        energy
        amplitude
        stability
        bifurcations

        validation_modes

        h_amplitude
        corrected_low_modal_amplitude
        low_modal_amplitude
        h_energy
        h_modal_energy_fraction

        periodicity_error
    end
    methods
        function obj = Dynamic_Data(Model,Continuation_Opts)
            if nargin == 1
                Continuation_Opts = struct([]);
            end
            obj = obj.update_continuation_opts(Continuation_Opts);
            obj.Dynamic_Model = Model;

            obj.num_solutions = 0;
            obj.solution_types = {};
        end
        %-----------------------------------------------------------------%
        function obj = update_continuation_opts(obj,Continuation_Opts)
            Default_Continuation_Opts = read_default_options("continuation");         
            obj.Continuation_Options = update_options(Default_Continuation_Opts,obj.Continuation_Options,Continuation_Opts);
        end 
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        % Creation
        %-----------------------------------------------------------------%
        function obj = add_backbone(obj,mode_num,type)
            if nargin == 2
                type = "rom";
            end
            Rom = obj.Dynamic_Model;
            solution_num = obj.num_solutions + 1;
            
            
            [t0,z0] = get_linear_solution(Rom,mode_num,type);
            
            

            coco_backbone(t0,z0,Rom,type,obj.Continuation_Options,solution_num);
            obj = obj.analyse_solution(solution_num);
            
            solution_type = "backbone";
            if type == "fom"
                solution_type = "fom: " + solution_type;
            end
            obj.solution_types{1,solution_num} = solution_type;
            obj.num_solutions = solution_num;

            obj.save_solution("coco",solution_num)
        end
        %-----------------------------------------------------------------%
        function obj = restart_backbone(obj,solution_num,orbit_num,type)
            next_solution_num = obj.num_solutions + 1;
            Rom = obj.Dynamic_Model;

            solution_type = obj.solution_types{1,solution_num};
            switch solution_type
                case "backbone"
                    backbone_type = "rom";
                case "fom: backbone"
                    backbone_type = "fom";
            end

            switch type
                case "po"
                    orbit = obj.get_orbit(solution_num,orbit_num);
                    t0 = orbit.tbp';
                    z0 = orbit.xbp';

                    coco_backbone(t0,z0,Rom,backbone_type,obj.Continuation_Options,next_solution_num);
            end

            obj = obj.analyse_solution(next_solution_num);

            obj.solution_types{1,next_solution_num} = solution_type;
            obj.num_solutions = next_solution_num;

            obj.save_solution("coco",next_solution_num)
        end
        %-----------------------------------------------------------------%
        function obj = add_full_order_backbone(obj,mode_num)
            Model = obj.Dynamic_Model.Model;
            if Model.system_type ~= "direct"
                error("Full order backbone curves are not supported for FE systems")
            end
            obj = obj.add_backbone(mode_num,"fom");
        end
        %-----------------------------------------------------------------%
        function obj = add_modal_frf(obj,Force_Data,Damping_Data)
            Rom = obj.Dynamic_Model;
            solution_num = obj.num_solutions + 1;
            switch Damping_Data.damping_type
                case "rayleigh"
                    damping = get_rayleigh_damping_matrix(Damping_Data,obj.Dynamic_Model.Model);
                    Nonconservative_Input.damping = damping;
                    Nonconservative_Input.amplitude = Force_Data.amplitude;
            end
            Eom_Input = Rom.get_solver_inputs("coco_backbone");

            [t0,z0] = get_forced_response(obj.Dynamic_Model,Eom_Input,Nonconservative_Input);

            coco_frf(t0,z0,Eom_Input,Nonconservative_Input,obj.Continuation_Options,solution_num);

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Validation
        %-----------------------------------------------------------------%
        function  obj = validate_solution(obj,solution_num,L_modes)
             H_ALGORITHM = "time_solution";
            
            if isstring(solution_num)
                if solution_num == "last"
                    solution_num = obj.num_solutions;
                end
            end
            data_path = obj.Dynamic_Model.data_path;
            load(data_path + "Static_Data.mat","Static_Data")
            Static_Data = Static_Data.add_validation_data(L_modes);
            save(data_path + "Static_Data.mat","Static_Data","-v7.3")

            L_modes = Static_Data.Dynamic_Validation_Data.current_L_modes;
            
            obj.validation_modes{1,solution_num} = L_modes;
            Validation_Rom = Reduced_System(Static_Data);

          
            switch H_ALGORITHM
                case "harmonic_balance"
                    obj = h_harmonic_balance(obj,Validation_Rom,solution_num);
                case "time_solution"
                    obj = h_time_solution(obj,Validation_Rom,solution_num);
            end
             obj.save_solution("h_prediciton");
        end
        %-----------------------------------------------------------------%
        function obj = get_periodicity_error(obj,solution_num,orbit_num)
            orbit = get_orbit(obj,solution_num,orbit_num);
            Rom = obj.Dynamic_Model;
            solution_periodicity_error(orbit,Rom)
        end

        %-----------------------------------------------------------------%
        % Analysis 
        %-----------------------------------------------------------------%
        function obj = analyse_solution(obj,solution_num)
            solution_name = "temp\dynamic_sol_" + solution_num; 
            bd = coco_bd_read(solution_name);

            solution_data  = coco_bd_col(bd, {"po.period", "LAB", "TYPE", "eigs","ENERGY","MAX(x)","MIN(x)"});
            
            %------------------------%
            % Labels
            orbit_labels = solution_data{1,2};
            
            %------------------------%
            % Frequency
            period = solution_data{1,1};
            orbit_frequency = 2*pi./period;
            
            %------------------------%
            % Energy
            orbit_energy = solution_data{1,5};
            
            %------------------------%
            % Bifurcations
            orbit_type = solution_data{1,3};
            orbit_bifurcation.bp_index = find(orbit_type == "BP");
            orbit_bifurcation.pd_index = find(orbit_type == "PD");
            orbit_bifurcation.ns_index = find(orbit_type == "TR");
            orbit_bifurcation.sn_index = find(orbit_type == "SN");
           
            %------------------------%
            % Stability
            eigs = solution_data{1,4};
            orbit_stability = max(abs(eigs),[],1);
            
            %------------------------%
            % Amplitude
             % sol = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(iSol));
            num_modes = size(solution_data{1,6},1)/2;
            min_disp = solution_data{1,6}(1:num_modes,:);
            max_disp = solution_data{1,7}(1:num_modes,:);
            orbit_amplitude = abs(max_disp - min_disp)/2;
            
 
            %------------------------%
            obj.solution_labels{1,solution_num} = orbit_labels;
            obj.frequency{1,solution_num} = orbit_frequency;
            obj.energy{1,solution_num} = orbit_energy;
            obj.amplitude{1,solution_num} = orbit_amplitude;
            obj.stability{1,solution_num} = orbit_stability;
            obj.bifurcations{1,solution_num} = orbit_bifurcation;

            obj.validation_modes{1,solution_num} = [];

            obj = obj.pre_allocate_validation_data(solution_num,0);
           
        end
        %-----------------------------------------------------------------%
        function obj = analyse_h_solution(obj,r,r_dot,h,h_dot,Validation_Analysis_Inputs,solution_num,orbit_num)
            
            
            get_amplitude = @(x) (max(x,[],2) - min(x,[],2))/2;
            %--
            h_amp = get_amplitude(h);
            
            %--
            Theta_Poly = Validation_Analysis_Inputs.Disp_Data.Theta_Poly;
            theta = Theta_Poly.evaluate_polynomial(r);

            L_disp_transform = Validation_Analysis_Inputs.L_disp_transform;
            g_l = L_disp_transform*theta;

            g_l = [r;g_l];
            g_amp = get_amplitude(g_l);

            %--
            q_l = h+g_l;
            q_amp = get_amplitude(q_l);

            %--
            num_points = size(r,2);
            num_h_modes = size(h,1);
            h_potential = zeros(num_h_modes,num_points);
            for iPoint = 1:num_points
                 h_stiffness = Validation_Analysis_Inputs.H_Stiffness_Poly.evaluate_polynomial(r(:,iPoint));
                 h_force = h_stiffness*h(:,iPoint);
                 h_potential(:,iPoint) = 1/2*(h_force.*h(:,iPoint));
            end
            
            potential_tilde = Validation_Analysis_Inputs.Potential.evaluate_polynomial(r);
            potential_hat = potential_tilde + sum(h_potential);
            %--

            %--
            [ke_mode_tilde,ke_condensed_tilde,ke_mode_hat,ke_condensed_hat] = h_kinetic_energy(r,r_dot,h,h_dot,Validation_Analysis_Inputs);
            %--
            % energy_test = obj.energy{1,solution_num}(orbit_num);
            % energy_tilde = potential_tilde + sum(ke_mode_tilde) + ke_condensed_tilde;
            energy_hat =  potential_hat + sum(ke_mode_hat) + ke_condensed_hat;
            % energy_fraction = (mean(energy_hat)-mean(energy_tilde))/mean(energy_hat);
            %--
            obj.h_amplitude{1,solution_num}(:,orbit_num) = h_amp;
            obj.corrected_low_modal_amplitude{1,solution_num}(:,orbit_num) = g_amp;
            obj.low_modal_amplitude{1,solution_num}(:,orbit_num) = q_amp;
            obj.h_energy{1,solution_num}(1,orbit_num) = mean(energy_hat);
            % obj.h_modal_energy_fraction{1,solution_num}(:,orbit_num) = max(ke_mode_fraction,[],2);
        end
        %-----------------------------------------------------------------%
        function orbit = get_orbit(obj,solution_num,orbit_num)
            data_path = obj.Dynamic_Model.data_path;
            orbit_labels = obj.solution_labels{1,solution_num};

            solution_name = data_path + "dynamic_sol_" + solution_num; 
            orbit = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(orbit_num));
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Helpers
        %-----------------------------------------------------------------%
        function save_solution(Dyn_Data,type,solution_num)
            data_path = Dyn_Data.Dynamic_Model.data_path;
            switch type
                case "coco"
                    solution_name = "dynamic_sol_" + solution_num + "\";
                    movefile("data\temp\" + solution_name,data_path + solution_name)
                case "h_prediction"

                case "update"

            end
            save(data_path + "Dyn_Data","Dyn_Data")
        end
        %-----------------------------------------------------------------%
        function obj = remove_solution(obj,solution_num)
            if isstring(solution_num)
                if solution_num == "all"
                    for iSol = 1:obj.num_solutions
                        obj = obj.remove_solution(1);
                    end
                    return
                end
            end
            %------
            solution_num = sort(solution_num,"descend");
            num_removed_solutions = length(solution_num);

            for iSol = 1:num_removed_solutions

                data_path = obj.Dynamic_Model.data_path;
                sol_path = data_path + "dynamic_sol_" + solution_num(iSol);
                if ~isfolder(sol_path)
                    warning("Solution cannot be removed as it does not exist")
                    return
                end
                rmdir(sol_path,'s')

                %------
                num_sols = obj.num_solutions;
                old_sol_names = (solution_num(iSol)+1) : num_sols;
                num_old_sols = length(old_sol_names);
                if ~isempty(old_sol_names)
                    new_sol_names = old_sol_names - 1;
                    for jSol = 1:num_old_sols
                        old_sol_path = data_path + "dynamic_sol_" + old_sol_names(jSol);
                        new_sol_path = data_path + "dynamic_sol_" + new_sol_names(jSol);
                        movefile(old_sol_path,new_sol_path)
                    end
                end

                %------
                dyn_data_properties = properties(obj);
                num_properties = length(dyn_data_properties);
                for iProperty = 1:num_properties
                    property = dyn_data_properties{iProperty};
                    if ~iscell(obj.(property))
                        continue
                    end

                    obj.(property)(:,solution_num(iSol)) = [];
                end
                obj.num_solutions = num_sols - 1;
            end
            save_solution(obj,"update")
        end
        %-----------------------------------------------------------------%
        function obj = pre_allocate_validation_data(obj,solution_num,num_orbits)
            num_r_modes = length(obj.Dynamic_Model.Model.reduced_modes);
            L_modes = obj.validation_modes{1,solution_num};
            num_h_modes = length(L_modes) + num_r_modes;
            obj.h_amplitude{1,solution_num} = zeros(num_h_modes,num_orbits);
            obj.corrected_low_modal_amplitude{1,solution_num} = zeros(num_h_modes,num_orbits);
            obj.low_modal_amplitude{1,solution_num} = zeros(num_h_modes,num_orbits);
            obj.h_energy{1,solution_num} = zeros(1,num_orbits);
            obj.h_modal_energy_fraction{1,solution_num} = zeros(num_h_modes,num_orbits);
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Overloading 
        %-----------------------------------------------------------------%
        function sz = size(obj)
            sz  = [obj.num_solutions,1];
        end
    end
        
end