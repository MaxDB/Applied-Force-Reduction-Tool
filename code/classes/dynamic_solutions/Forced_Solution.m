classdef Forced_Solution < Dynamic_Solution
    properties
        frequency
        energy
        amplitude
        stability
        bifurcations

        Force_Data
        Damping_Data

        additional_dynamic_output
    end
    
    methods
        function obj = Forced_Solution(Rom,FRF_Settings)
            obj@Dynamic_Solution(FRF_Settings.Continuation_Opts)
            
            Force_Data = FRF_Settings.Force_Data;
            Damping_Data = FRF_Settings.Damping_Data;
            type = FRF_Settings.type;
            Add_Ouput = FRF_Settings.Additional_Output;
            solution_num = FRF_Settings.solution_num;
            initial_conditions = FRF_Settings.initial_condition;
            
            obj.Force_Data = Force_Data;
            obj.Damping_Data = Damping_Data;

            continuation_variable = Force_Data.continuation_variable;
            Nonconservative_Input = obj.get_nonconservative_input(Rom);
            if isfield(FRF_Settings,"z0")
                Nonconservative_Input.z0 = FRF_Settings.z0;
            end


            switch continuation_variable
                case "amplitude"
                    [t0,z0,F0] = get_forced_linear_solution(Rom,Nonconservative_Input,continuation_variable);
                    p0 = F0;

                    Sol_Type.frequency = Force_Data.frequency;
                case "frequency"
                    if isempty(initial_conditions)
                        % period_range = obj.Continuation_Options.parameter_range;
                        % p0 = period_range(2);
                        p0 = Force_Data.frequency;
                        T0 = 2*pi/Force_Data.frequency;
                        [t0,z0] = get_forced_response(Rom,Nonconservative_Input,T0,type);
                    else
                        [t0,z0,p0] = initial_conditions{:};
                    end

                    Sol_Type.amplitude = Force_Data.amplitude;
            end


            coco_forced_response(t0,z0,p0,Rom,type,obj.Continuation_Options,solution_num,Add_Ouput,Nonconservative_Input);
            
            
            obj = obj.analyse_solution(solution_num,Add_Ouput);
            
            Sol_Type.orbit_type = "forced";
            Sol_Type.model_type = type;
     
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
            % F_Data.shape = Applied_Force_Data.shape;
            % F_Data.type = "shape";
            
            switch Damp_Data.damping_type
                case "rayleigh"
                    damping = get_rayleigh_damping_matrix(Damp_Data,Model);
                    Nonconservative_Input.damping = damping;
                    Nonconservative_Input.damping_type = "matrix";
                case "modal"
                    damping_factors = Damp_Data.modal_damping_factors;
                    damping = diag(damping_factors);
                    Nonconservative_Input.damping = damping;
                    Nonconservative_Input.damping_type = "matrix";
                case "physical"
                    Nonconservative_Input.damping = Damp_Data.damping_matrix;
                    Nonconservative_Input.damping_type = "matrix";
                case "nonlinear_rayleigh"
                    Nonconservative_Input.damping_factors = [Damp_Data.mass_factor,Damp_Data.stiffness_factor];
                    Nonconservative_Input.damping_type = "nonlinear_rayleigh";
                otherwise
                    error("Damping model unsupported / not recognised ")
            end
            
            continuation_variable = F_Data.continuation_variable;

            switch F_Data.type
                case "modal"
                    mode_map = F_Data.mode_number == Model.reduced_modes;
                    
                    shape = Model.mass*Model.reduced_eigenvectors;
                    shape = shape(:,mode_map);
                    Nonconservative_Input.amplitude_shape = shape;
                case "point"
                    num_dofs = Model.num_dof + length(Model.dof_boundary_conditions);
                    dof_map = zeros(num_dofs,1);
                    dof_map(F_Data.dof) = 1;
                    dof_map(Model.dof_boundary_conditions) = [];

                    Nonconservative_Input.amplitude_shape = dof_map;
                case "shape"
                    Nonconservative_Input.amplitude_shape = F_Data.shape;
                case "uniform"
                    num_dof = Rom.Model.num_dof +  length(Rom.Model.dof_boundary_conditions);
                    num_dimensions = get_num_node_dimensions(Rom.Model);
                    num_nodes = num_dof/num_dimensions;
                    shape = zeros(num_dof,1);
                    direction_index = (0:(num_nodes-1))*num_dimensions + F_Data.direction;
                    shape(direction_index) = 1;
                    shape(Rom.Model.dof_boundary_conditions) = [];
                    shape = shape/norm(shape);

                    Nonconservative_Input.amplitude_shape = shape;

            end
            Nonconservative_Input.force_type = "shape";
            Nonconservative_Input.continuation_variable = continuation_variable;

            switch continuation_variable
                case "amplitude"
                    Nonconservative_Input.frequency = F_Data.frequency;
                    Nonconservative_Input.force_points = F_Data.force_points;
                case "frequency"
                    Nonconservative_Input.amplitude = F_Data.amplitude;
            end

        end
        %-----------------------------------------------------------------%
        function obj = analyse_solution(obj,solution_num,Add_Output)
            Analysis_Output = analyse_solution@Dynamic_Solution(obj,solution_num,Add_Output);
            obj = Dynamic_Solution.update_dynamic_solution(obj,Analysis_Output);
        end
    end

end