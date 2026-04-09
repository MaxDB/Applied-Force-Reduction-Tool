classdef FRF_To_BB_Solution < Dynamic_Solution
    properties
        frequency
        energy
        amplitude
        stability
        bifurcations

        Force_Data
        Damping_Data

        epsilon

        additional_dynamic_output
    end
    
    methods
        function obj = FRF_To_BB_Solution(Rom,FRF_Settings)
            obj@Dynamic_Solution(FRF_Settings.Continuation_Opts)
            
            Force_Data = FRF_Settings.Force_Data;
            Damping_Data = FRF_Settings.Damping_Data;
            type = FRF_Settings.type;
            Add_Ouput = FRF_Settings.Additional_Output;
            solution_num = FRF_Settings.solution_num;
            initial_conditions = FRF_Settings.initial_condition;
            
            obj.Force_Data = Force_Data;
            obj.Damping_Data = Damping_Data;

            Nonconservative_Input = obj.get_nonconservative_input(Rom);
            if isfield(FRF_Settings,"z0")
                Nonconservative_Input.z0 = FRF_Settings.z0;
            end
            
            [t0,z0,p0] = initial_conditions{:};
            period = t0(end);
            Sol_Type.frequency = 2*pi/period;
            Sol_Type.amplitude = Force_Data.amplitude;

            if p0 == 0
                T0 = 2*pi/Force_Data.frequency;
                error("Need to take into account phase difference between applied force and response")
                [t0,z0] = get_forced_response(Rom,Nonconservative_Input,T0,type);
            end

            coco_frf_to_bb(t0,z0,p0,period,Force_Data.amplitude,Rom,type,obj.Continuation_Options,solution_num,Add_Ouput,Nonconservative_Input);
            
            
            obj = obj.analyse_solution(solution_num,Add_Ouput);
            
            Sol_Type.orbit_type = "frf_to_bb";
            Sol_Type.model_type = type;
     
            obj.Solution_Type = Sol_Type;
        end
        %-----------------------------------------------------------------%
        function Nonconservative_Input = get_nonconservative_input(obj,Rom)
            F_Data = obj.Force_Data;
            Damp_Data = obj.Damping_Data;
            Model = Rom.Model;
            
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
                    dof_map(Model.node_mapping(:,1) == F_Data.dof) = 1;
                    Nonconservative_Input.amplitude_shape = dof_map;
                case "shape"
                    Nonconservative_Input.amplitude_shape = F_Data.shape;
                    
            end

            Nonconservative_Input.amplitude = F_Data.amplitude;
            Nonconservative_Input.frequency = F_Data.frequency;
            Nonconservative_Input.force_type = F_Data.type;
            Nonconservative_Input.continuation_variable = continuation_variable;

        end
        %-----------------------------------------------------------------%
        function obj = analyse_solution(obj,solution_num,Add_Output)
            Analysis_Output = analyse_solution@Dynamic_Solution(obj,solution_num,Add_Output,"epsilon");
            obj = Dynamic_Solution.update_dynamic_solution(obj,Analysis_Output);
            obj.epsilon = Analysis_Output.extra;
        end
    end

end