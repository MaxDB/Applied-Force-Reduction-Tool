classdef Sparse_Stiffness_Pointer
    properties
        degrees_of_freedom
        number_of_loadcases
        
        stiffness_name
        file_type
        data_dir
        object_name
    end

    methods
        function obj = Sparse_Stiffness_Pointer(dofs,varargin)
            obj_name = "Sparse_Stiffness";
            

            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            stiffness_path = "temp\stiffness";
            job_id = [];

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "job"
                        job_id = keyword_values{arg_counter};
                    case "path"
                        stiffness_path = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------%



            
            obj.object_name = obj_name;
            
            if ~isempty(job_id)
                stiffness_path = stiffness_path + "\" + job_id;
            end

            obj.data_dir = stiffness_path;

            if isfolder(stiffness_path)
                obj = obj.load();
                if obj.degrees_of_freedom ~= dofs
                    error("Incompatable stiffness tensors")
                end
                return
            else
                mkdir(stiffness_path)
            end

            obj.degrees_of_freedom = dofs;
            obj.number_of_loadcases = 0;

            
            obj.stiffness_name = "STIF";
            obj.file_type = ".mtx";
            
            
            
            obj.save;
            
        end
        %------
        function obj = add_loadcases(obj,base_file_name,step_list)
            num_loadcases = obj.number_of_loadcases;
            num_added_loadcases = size(step_list,2);
            for iLoad = 1:num_added_loadcases
                old_file_name = base_file_name + step_list(iLoad) + obj.file_type;
                new_file_name = obj.data_dir + "\" + obj.stiffness_name  + (num_loadcases + iLoad) + obj.file_type;
                movefile(old_file_name,new_file_name)
            end
            
            obj.number_of_loadcases = num_loadcases + num_added_loadcases;
            obj.save;
        end
        %------
        function obj = combine_data(obj,obj_two)
            if obj.degrees_of_freedom ~= obj_two.degrees_of_freedom
                error("Incompatable stiffness tensors")
            end
            if obj.data_dir == obj_two.data_dir
                error("Cannot add a stiffness point to itself")
            end
            
            destination_path = obj.data_dir;
            current_path = obj_two.data_dir;
            num_old_loadcases = obj.number_of_loadcases;
            num_added_loadcases = obj_two.number_of_loadcases;
            
            for iLoad = 1:num_added_loadcases
                new_load_label = iLoad + num_old_loadcases;
                stiffness_file = current_path + "\" + obj_two.stiffness_name + iLoad + obj_two.file_type;
                new_stiffness_file = destination_path + "\" + obj.stiffness_name + new_load_label + obj.file_type;
                movefile(stiffness_file,new_stiffness_file)
            end
    

            rmdir(current_path,'s');
            obj.number_of_loadcases = num_old_loadcases + num_added_loadcases;
        end
        %------
        function save(Sparse_Stiffness)
            file_dir = Sparse_Stiffness.data_dir + "\" +Sparse_Stiffness.object_name+"_data";
            save(file_dir,Sparse_Stiffness.object_name);
        end
        %------
        function Sparse_Stiffness = load(obj) %#ok<STOUT>
            load(obj.data_dir + "\" + obj.object_name+ "_data.mat",obj.object_name);
        end
        %------



        %-----------------------
        % overloading
        %------
        function obj = cat(dim,varargin)
            if dim ~= 3
                error("Not implemented")
            end

            num_arrays = size(varargin,2);
            obj = varargin{1};
            
            for iArray = 2:num_arrays
                obj = obj.combine_data(varargin{iArray});
            end
            obj.save;
        end
        %------
        function obj_size = size(obj,dim)
            obj_size = [obj.degrees_of_freedom,obj.degrees_of_freedom,obj.number_of_loadcases];
            if nargin == 2
               obj_size = obj_size(dim); 
            end
        end
    end


    end
   