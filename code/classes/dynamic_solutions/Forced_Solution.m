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
            Nonconservative_Input = obj.get_nonconservative_input(Rom.Model);
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
        function Nonconservative_Input = get_nonconservative_input(obj,Model)
            F_Data = obj.Force_Data;
            Damp_Data = obj.Damping_Data;
            
            switch Damp_Data.damping_type
                case "rayleigh"
                    damping = get_rayleigh_damping_matrix(Damp_Data,Model);
                    Nonconservative_Input.damping = damping;
            end
            continuation_variable = F_Data.continuation_variable;

            switch F_Data.type
                case "modal"
                    mode_map = F_Data.mode_number == Model.reduced_modes;
                    Nonconservative_Input.mode_map = mode_map;
                case "point force"
                    num_dofs = Model.num_dof;
                    dof_map = zeros(num_dofs,1);
                    if isempty(Model.dof_boundary_conditions)
                        dof_map(F_Data.dof) = 1;
                    else
                        dof_map(Model.node_mapping(:,1) == F_Data.dof) = 1;
                    end

                    Nonconservative_Input.amplitude_shape = dof_map;
            end
            Nonconservative_Input.force_type = F_Data.type;
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