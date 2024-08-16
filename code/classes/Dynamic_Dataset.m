classdef Dynamic_Dataset
        %Creates and stores dynamic results: periodic orbits, backbone curves, FRFs ect... 
    properties
        num_solutions
        solution_types
        
        Dynamic_Model
        Additional_Output
    end
    
    methods
        function obj = Dynamic_Dataset(Model)
            obj.Dynamic_Model = Model;
            
            obj.Additional_Output.type = "none";
            obj.num_solutions = 0;
            obj.solution_types = {};
        end
        %-----------------------------------------------------------------%
        function obj = add_additional_output(obj,Additional_Output)
            obj.Additional_Output = Additional_Output;
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Creation
        %-----------------------------------------------------------------%
        function obj = add_backbone(obj,mode_num,varargin)
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Continuation_Opts = struct([]);
            type = "rom";
            initial_condition = [];

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "type"
                        type = keyword_values{arg_counter};
                    case "opts"
                        Continuation_Opts = keyword_values{arg_counter};
                    case "ic"
                        initial_condition = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------%
            BB_Settings.mode_num = mode_num;
            BB_Settings.type = type;
            BB_Settings.solution_num = obj.num_solutions + 1;
            BB_Settings.Continuation_Opts = Continuation_Opts;
            BB_Settings.Additional_Output = obj.Additional_Output;
            BB_Settings.initial_condition = initial_condition;

            Rom = obj.Dynamic_Model;
            BB_Sol = Backbone_Solution(Rom,BB_Settings);

            obj.num_solutions = obj.num_solutions + 1;
            obj.solution_types{obj.num_solutions} = BB_Sol.Solution_Type;
            obj.solution_types{obj.num_solutions}.validated = false;
            obj.save_solution(BB_Sol)
        end
        %-----------------------------------------------------------------%        
        function obj = restart_point(obj,solution_num,orbit_num,type,varargin)
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Continuation_Opts = struct([]);

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "opts"
                        Continuation_Opts= keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------%
            Rom = obj.Dynamic_Model;
            
            switch type
                case "IC"
                    point_index = obj.get_special_point(solution_num,type);
            end

            if isstring(orbit_num)
                if orbit_num == "all"
                    orbit_num = 1:size(point_index,1);
                end
            end
            num_orbits = length(orbit_num);
            for iOrbit = 1:num_orbits
                next_solution_num = obj.num_solutions + 1;

                Solution_Type = obj.solution_types{1,solution_num};
                model_type = Solution_Type.model_type;

                switch type
                    case "po"
                        orbit = obj.get_orbit(solution_num,orbit_num(iOrbit));
                        t0 = orbit.tbp';
                        z0 = orbit.xbp';

                        obj = obj.add_backbone(0,"opts",Continuation_Opts,"ic",{t0,z0});
                    case "bp"

                    case "IC"
                        orbit = obj.get_orbit(solution_num,point_index(orbit_num(iOrbit)));

                        t0 = orbit.tbp';
                        z0 = orbit.xbp';
                        p0 = Solution_Type.frequency;

                        FRF_Sol = obj.load_solution(solution_num);
                        Force_Data = FRF_Sol.Force_Data;
                        Damping_Data = FRF_Sol.Damping_Data;

                        Force_Data.continuation_variable = "frequency";
                        Force_Data = rmfield(Force_Data,"frequency");
                        Force_Data.amplitude = orbit.p;

                        obj = obj.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"ic",{t0,z0,p0});
                        % coco_forced_response(t0,z0,p0,Rom,model_type,obj.Continuation_Options,next_solution_num,obj.Additional_Output,Nonconservative_Inputs);
                        % 
                        % obj = obj.analyse_solution(next_solution_num);
                        % 
                        % Solution_Type = rmfield(Solution_Type,"frequency");
                        % Solution_Type.amplitude = orbit.p;
                        % obj.solution_types{1,next_solution_num} = Solution_Type;
                        % obj.num_solutions = next_solution_num;
                        % 
                        % obj.save_solution("coco_frf",next_solution_num,Nonconservative_Inputs)
                    case "force"
                        orbit = obj.get_orbit(solution_num,orbit_num(iOrbit));

                        t0 = orbit.tbp';
                        z0 = orbit.xbp';
                        p0 = Solution_Type.amplitude;



                        data_path = obj.Dynamic_Model.data_path;
                        solution_name = data_path + "dynamic_sol_" + solution_num;
                        load(solution_name + "\Nonconservative_Inputs.mat","Nonconservative_Inputs")

                        Nonconservative_Inputs.continuation_variable = "amplitude";
                        Nonconservative_Inputs.frequency = obj.frequency{1,solution_num}(orbit_num(iOrbit));

                        Nonconservative_Inputs = rmfield(Nonconservative_Inputs,"amplitude");

                        coco_forced_response(t0,z0,p0,Rom,model_type,obj.Continuation_Options,next_solution_num,obj.Additional_Output,Nonconservative_Inputs);
                        obj = obj.analyse_solution(next_solution_num);

                        Solution_Type = rmfield(Solution_Type,"amplitude");
                        Solution_Type.frequency = obj.frequency{1,solution_num}(orbit_num(iOrbit));
                        obj.solution_types{1,next_solution_num} = Solution_Type;
                        obj.num_solutions = next_solution_num;

                        obj.save_solution("coco_frf",next_solution_num,Nonconservative_Inputs)
                end


            end
        end
        %-----------------------------------------------------------------%
        function obj = add_forced_response(obj,Force_Data,Damping_Data,varargin)
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Continuation_Opts = struct([]);
            type = "rom";
            initial_condition = [];

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "type"
                        type = keyword_values{arg_counter};
                    case "opts"
                        Continuation_Opts = keyword_values{arg_counter};
                    case "ic"
                        initial_condition = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------%


            % obj.num_solutions = obj.num_solutions + 1;
            % obj.solution_types{obj.num_solutions} = BB_Sol.Solution_Type;
            % obj.solution_types{obj.num_solutions}.validated = false;
            % obj.save_solution(BB_Sol)

            FRF_Settings.Force_Data = Force_Data;
            FRF_Settings.Damping_Data = Damping_Data;
            FRF_Settings.Continuation_Opts = Continuation_Opts;
            FRF_Settings.solution_num = obj.num_solutions + 1;
            FRF_Settings.Additional_Output = obj.Additional_Output;
            FRF_Settings.type = type;
            FRF_Settings.initial_condition = initial_condition;

            Rom = obj.Dynamic_Model;
            FRF_Sol = Forced_Solution(Rom,FRF_Settings);

           
            obj.num_solutions = obj.num_solutions + 1;
            obj.solution_types{obj.num_solutions} = FRF_Sol.Solution_Type;
            obj.solution_types{obj.num_solutions}.validated = false;
            obj.save_solution(FRF_Sol)
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Validation
        %-----------------------------------------------------------------%
        function  obj = validate_solution(obj,solution_num,L_modes)
            if isstring(solution_num)
                if solution_num == "last"
                    solution_num = obj.num_solutions;
                end
            end

            Rom = obj.Dynamic_Model;

            

            BB_Sol = obj.load_solution(solution_num);

            if isstring(L_modes) && L_modes == "all"
                switch Rom.Model.system_type
                    case "indirect"
                        L_modes = Rom.Model.low_frequency_modes;
                    case "direct"
                        L_modes = 1:Rom.Model.num_dof;
                end
            end
            Validated_BB_Settings.solution_num = solution_num;
            Validated_BB_Settings.L_modes = L_modes;

            Validated_BB_Sol = Validated_Backbone_Solution(Rom,BB_Sol,Validated_BB_Settings);

            obj.solution_types{obj.num_solutions}.validated = true;
            obj.save_solution(Validated_BB_Sol)
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Helpers
        %-----------------------------------------------------------------%
        function update_dyn_data(Dyn_Data)
            data_path = Dyn_Data.Dynamic_Model.data_path;
            save(data_path + "Dyn_Data","Dyn_Data")
        end
        %-----------------------------------------------------------------%
        function save_solution(Dyn_Data,Solution)
            data_path = Dyn_Data.Dynamic_Model.data_path;
            solution_num = Dyn_Data.num_solutions;

            solution_name = "dynamic_sol_" + solution_num + "\";
            solution_path = data_path + solution_name;

            switch class(Solution)
                case {"Backbone_Solution","Forced_Solution"}
                    file_name = "Sol_Data";
                    movefile("data\temp\" + solution_name,solution_path);
                case "Validated_Backbone_Solution"
                    file_name = "Sol_Data_Validated";
            end 

            save(solution_path + file_name,"Solution")
            Dyn_Data.update_dyn_data;
        end
        %-----------------------------------------------------------------%
        function Solution = load_solution(obj,solution_num,type)
            if nargin == 2
                type = "Sol_Data";
            end

            switch type
                case "validation"
                    type = "Sol_Data_Validated";
            end
            data_path = obj.Dynamic_Model.data_path;

            solution_name = "dynamic_sol_" + solution_num + "\";
            solution_path = data_path + solution_name;
            load(solution_path + type,"Solution")
            
            switch type
                case "Sol_Data"
                     Solution.Solution_Type = obj.solution_types{solution_num};
                case "Sol_Data_Validated"
           
            end
            
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
                else
                    rmdir(sol_path,'s')
                end
                
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
                obj.num_solutions = num_sols - 1;
                obj.solution_types(:,solution_num(iSol)) = [];
            end
            obj.update_dyn_data;
        end
        %-----------------------------------------------------------------%
        function [orbit,validation_orbit] = get_orbit(obj,solution_num,orbit_num,validated)
            if nargin == 3
                validated = 0;
            end
            Solution = obj.load_solution(solution_num);
            
            data_path = obj.Dynamic_Model.data_path;
            orbit_labels = Solution.orbit_labels;

            solution_name = data_path + "dynamic_sol_" + solution_num; 
            
            
            if isscalar(orbit_num)
                orbit = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(orbit_num));
            else
                num_orbits = length(orbit_num);
                orbit = cell(num_orbits,1);
                for iOrbit = 1:num_orbits
                    orbit{iOrbit,1} = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(orbit_num(iOrbit)));
                end
            end

            if validated
                load(solution_name + "\sol" + orbit_num + "_v.mat","validation_orbit");
            else
                validation_orbit = [];
            end
        end
        %-----------------------------------------------------------------%
        function point_index = get_special_point(obj,solution_num,point_type)
            data_path = obj.Dynamic_Model.data_path;
            solution_name = data_path + "dynamic_sol_" + solution_num;
            bd = coco_bd_read(solution_name);
            orbit_type = coco_bd_col(bd, "TYPE");
            point_index = find(orbit_type == point_type);
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Overloading 
        %-----------------------------------------------------------------%
        function sz = size(obj,index)
            sz  = [obj.num_solutions,1];
            if nargin == 2
                sz = sz(index);
            end
        end
    end

end