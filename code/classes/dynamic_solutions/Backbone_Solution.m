classdef Backbone_Solution < Dynamic_Solution
    properties
        frequency
        energy
        amplitude
        stability
        bifurcations

        additional_dynamic_output
    end
    
    methods
        function obj = Backbone_Solution(Rom,BB_Settings)
            obj@Dynamic_Solution(BB_Settings.Continuation_Opts)
            

            mode_num = BB_Settings.mode_num;
            type = BB_Settings.type;
            [t0,z0] = get_linear_solution(Rom,mode_num,type);
            
            solution_num = BB_Settings.solution_num;
            Add_Ouput = BB_Settings.Additional_Output;
            coco_backbone(t0,z0,Rom,type,obj.Continuation_Options,solution_num,Add_Ouput);
            obj = obj.analyse_solution(solution_num,Add_Ouput);
            
            Sol_Type.orbit_type = "free";
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