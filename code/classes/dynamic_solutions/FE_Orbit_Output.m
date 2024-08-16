classdef FE_Orbit_Output
    properties
        orbit_labels
        solution_num
        
        num_orbits
        
        fe_output
        fe_output_type
    end

    methods
        function obj = FE_Orbit_Output(Dyn_Data,fe_output_type,solution_num,orbit_ids)
            obj.orbit_labels = orbit_ids;
            obj.num_orbits = size(orbit_ids,2);
            obj.solution_num = solution_num;
            obj.fe_output_type = fe_output_type;
            
            fe_data = obj.get_fe_output(Dyn_Data,orbit_ids);
           
            obj.fe_output = fe_data;

        end
        %-----------------------------------------------------------------%
        function obj = add_orbits(obj,Dyn_Data,orbit_ids)
            current_orbit_ids = obj.orbit_labels;
            new_orbit_ids = setdiff(orbit_ids,current_orbit_ids);
            if isempty(new_orbit_ids)
                return
            end
            fe_data = obj.get_fe_output(Dyn_Data,new_orbit_ids);
            obj.orbit_labels = [obj.orbit_labels,new_orbit_ids];
            obj.num_orbits = size(obj.orbit_labels,2);
            obj.fe_output = [obj.fe_output,fe_data];
        end
        %-----------------------------------------------------------------%
        function fe_data = get_fe_output(obj,Dyn_Data,orbit_ids)
            orbit = Dyn_Data.get_orbit(obj.solution_num,orbit_ids);
            Rom = Dyn_Data.Dynamic_Model;
            switch obj.fe_output_type
                case "periodicity"
                    fe_data = solution_periodicity_error(orbit,Rom);
                case "stress"
                    Add_Output = Dyn_Data.Additional_Output;
                    fe_data = solution_max_disp_stress(orbit,Rom,Add_Output);
            end
        end
        %-----------------------------------------------------------------%
    end

end