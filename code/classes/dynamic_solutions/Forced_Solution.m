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

            switch Damping_Data.damping_type
                case "rayleigh"
                    damping = get_rayleigh_damping_matrix(Damping_Data,Rom.Model);
                    Nonconservative_Input.damping = damping;
            end


            mode_map = Force_Data.mode_number==Rom.Model.reduced_modes;
            Nonconservative_Input.mode_map = mode_map;
            Nonconservative_Input.force_type = Force_Data.type;
            Nonconservative_Input.continuation_variable = continuation_variable;

            switch continuation_variable
                case "amplitude"
                    Nonconservative_Input.frequency = Force_Data.frequency;
                    Nonconservative_Input.force_points = Force_Data.force_points;


                    [t0,z0,F0] = get_forced_linear_solution(Rom,Nonconservative_Input,continuation_variable);
                    p0 = F0;

                    Sol_Type.frequency = Force_Data.frequency;
                case "frequency"
                    Nonconservative_Input.amplitude = Force_Data.amplitude;
                    if isempty(initial_conditions)

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
        function obj = analyse_solution(obj,solution_num,Add_Output)
            Analysis_Output = analyse_solution@Dynamic_Solution(obj,solution_num,Add_Output);
            obj = Dynamic_Solution.update_dynamic_solution(obj,Analysis_Output);
        end
    end

end