classdef Full_Order_Forced_Solution
    properties
        Solution_Type
        num_orbits

        frequency
        energy
        amplitude

        periodicity
        stability

        Force_Data
        Damping_Data

        additional_dynamic_output
    end
    
    methods
        function obj = Full_Order_Forced_Solution(Rom,FRF_Settings)
            
            Force_Data = FRF_Settings.Force_Data;
            Damping_Data = FRF_Settings.Damping_Data;
            Add_Ouput = FRF_Settings.Additional_Output;
            solution_num = FRF_Settings.solution_num;
            initial_conditions = FRF_Settings.initial_condition;
            
            obj.Force_Data = Force_Data;
            obj.Damping_Data = Damping_Data;

            Nonconservative_Input = obj.get_nonconservative_input(Rom);
            


            
            % get_full_order_forced_response(Rom.Model,Nonconservative_Input,solution_num);
            obj.frequency = Nonconservative_Input.frequency;
            obj.num_orbits = length(obj.frequency);

            obj = obj.analyse_solution(Rom.Model,solution_num,Add_Ouput);
            

            Sol_Type.orbit_type = "forced";
            Sol_Type.model_type = "fom";
     
            obj.Solution_Type = Sol_Type;
        end
        %-----------------------------------------------------------------%
        function Nonconservative_Input = get_nonconservative_input(obj,Rom)
            F_Data = obj.Force_Data;
            Damp_Data = obj.Damping_Data;
            Model = Rom.Model;

            Applied_Force_Data = Rom.Applied_Force_Data;
            if ~isempty(Applied_Force_Data)
                error("")
            end

            switch Damp_Data.damping_type
                case "rayleigh"
                    Nonconservative_Input.alpha = Damp_Data.mass_factor;
                    Nonconservative_Input.beta = Damp_Data.stiffness_factor;
            end
            

            switch F_Data.type
                case "modal"
                    mode_map = F_Data.mode_number == Model.reduced_modes;
                    
                    shape = Model.mass*Model.reduced_eigenvectors;
                    shape = shape(:,mode_map);
                    Nonconservative_Input.force_shape = shape;
                case "point"
                    num_dofs = Model.num_dof + length(Model.dof_boundary_conditions);
                    dof_map = zeros(num_dofs,1);
                    dof_map(F_Data.dof) = 1;
                    dof_map(Model.dof_boundary_conditions) = [];

                    Nonconservative_Input.force_shape = dof_map;
                case "shape"
                    Nonconservative_Input.force_shape = F_Data.shape;
                case "uniform"
                    num_dof = Rom.Model.num_dof +  length(Rom.Model.dof_boundary_conditions);
                    num_dimensions = get_num_node_dimensions(Rom.Model);
                    num_nodes = num_dof/num_dimensions;
                    shape = zeros(num_dof,1);
                    direction_index = (0:(num_nodes-1))*num_dimensions + F_Data.direction;
                    shape(direction_index) = 1;
                    shape(Rom.Model.dof_boundary_conditions) = [];
                    shape = shape/norm(shape);

                    Nonconservative_Input.force_shape = shape;

            end
            Nonconservative_Input.force_type = "shape";
            % all_dofs = Rom.Model.num_dof + numel(Rom.Model.dof_boundary_conditions);
            % node_map = 1:all_dofs;
            % node_map(Model.dof_boundary_conditions) = [];
            
            % shape = zeros(all_dofs,1);
            % shape(node_map) = Nonconservative_Input.force_shape;
            % Nonconservative_Input.force_shape = shape;
            Nonconservative_Input.amplitude = F_Data.amplitude;
            Nonconservative_Input.frequency = F_Data.frequency;

        end
        %-----------------------------------------------------------------%
        function obj = analyse_solution(obj,Model,solution_num,Add_Output)
            
            data_path = split(Model.get_data_path,"\");
            data_path = join(data_path(1:2),"\") + "\dynamic_sol_" + solution_num;

   
            
            

            orbit_energy = zeros(1,obj.num_orbits);
            orbit_periodicity = zeros(1,obj.num_orbits);
            orbit_stability = ones(1,obj.num_orbits);
            orbit_add_output = nan(1,obj.num_orbits);
            for iOrbit = 1:obj.num_orbits
                Orbit_Data = load(data_path+"\sol" + iOrbit);
                Orbit_Data = Orbit_Data.Sol_Data; %FIx

                %------------------------%
            

                %------------------------%
                % Energy
                orbit_energy(iOrbit) = max(Orbit_Data.energy);


                %------------------------%
                % Amplitude
                % sol = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(iSol));
                % num_modes = size(solution_data{1,6},1)/2;
                % min_disp = solution_data{1,6}(1:num_modes,:);
                % max_disp = solution_data{1,7}(1:num_modes,:);
                % orbit_amplitude = abs(max_disp - min_disp)/2;
                % 
                % switch Add_Output.output
                %     case "physical displacement"
                %         additional_output  = coco_bd_col(bd, "DISP");
                % end

                %------------------------%
                % Convergence
                orbit_periodicity(iOrbit) = Orbit_Data.periodicity_error;

                %------------------------%
                switch Add_Output.output
                    case "physical displacement"
                        orbit_add_output(iOrbit) = Add_Output.output_func(Orbit_Data.disp);
                end
            end
            obj.energy = orbit_energy;
            obj.periodicity = orbit_periodicity;
            obj.stability = orbit_stability;
            obj.additional_dynamic_output = orbit_add_output;
      
        end
    end

end