classdef Validated_Backbone_Solution < Validated_Solution
    properties

    end

    methods
        function obj = Validated_Backbone_Solution(Rom,BB_Sol,Validated_BB_Settings)
            obj@Validated_Solution(Rom,BB_Sol,Validated_BB_Settings)
        end
        
        %-----------------------------------------------------------------%
        function obj = analyse_h_solution(obj,Displacment,Velocity,Eom_Terms,orbit_stability,Validation_Analysis_Inputs,orbit_num)
            Analysis_Output = analyse_h_solution@Validated_Solution(obj,Displacment,Velocity,Eom_Terms,orbit_stability,Validation_Analysis_Inputs,orbit_num);
            obj = Validated_Solution.update_validated_solution(obj,Analysis_Output);
        end
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%

    end

end