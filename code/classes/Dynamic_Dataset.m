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
            %-------------------------------------------------------------------------%
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Continuation_Opts = struct([]);
            type = "rom";

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "type"
                        type = keyword_values{arg_counter};
                    case "opts"
                        Continuation_Opts= keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------%
            BB_Settings.mode_num = mode_num;
            BB_Settings.type = type;
            BB_Settings.solution_num = obj.num_solutions + 1;
            BB_Settings.Continuation_Opts = Continuation_Opts;
            BB_Settings.Additional_Output = obj.Additional_Output;
            
            Rom = obj.Dynamic_Model;
            BB_Sol = Backbone_Solution(Rom,BB_Settings);

            obj.num_solutions = obj.num_solutions + 1;
            obj.solution_types{obj.num_solutions} = BB_Sol.solution_type;
            obj.save_solution(BB_Sol)
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

            movefile("data\temp\" + solution_name,solution_path);

            save(solution_path + "Sol_Data","Solution")
            Dyn_Data.update_dyn_data;
        end
        %-----------------------------------------------------------------%
        function Solution = load_solution(obj,solution_num)
            data_path = obj.Dynamic_Model.data_path;

            solution_name = "dynamic_sol_" + solution_num + "\";
            solution_path = data_path + solution_name;
            load(solution_path + "Sol_Data","Solution")
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
                obj.solution_types(solution_num(iSol),:) = [];
            end
            obj.update_dyn_data;
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