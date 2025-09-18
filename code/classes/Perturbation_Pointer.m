classdef Perturbation_Pointer
    properties
        number_of_loadcases
        perturbed_modes
        degrees_of_freedom
        num_validation_modes

        perturbation_name
        file_type
        data_dir
        object_name
    end

    methods
        function obj = Perturbation_Pointer(Model,perturbation_disp,varargin)
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            perturbation_path = "temp\perturbation";

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "path"
                        perturbation_path = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------

            obj.file_type = ".mat";
            obj.object_name = "Perturbation_Pointer";
            obj.perturbation_name = "perturb";
            obj.data_dir = perturbation_path;

            if isfolder(obj.data_dir)
                rmdir(obj.data_dir,"s");
            end
            mkdir(obj.data_dir);
            
            if isempty(perturbation_disp)
                obj.number_of_loadcases = 0;
                obj.num_validation_modes = 0;
                obj.perturbed_modes = [];
            else
                obj.number_of_loadcases = size(perturbation_disp,3);
                obj.num_validation_modes = size(perturbation_disp,2);
                obj.perturbed_modes = [Model.reduced_modes,Model.low_frequency_modes];
            end
            obj.degrees_of_freedom = Model.num_dof;

         
            

            obj.save_data(perturbation_disp);
            obj.save();
        end
        %-----
        function perturbation_displacements = get_displacement(obj,h_modes)
            if ~all(ismember(h_modes,obj.perturbed_modes))
                perturbation_displacements = [];
                return
            end
            num_h_modes = size(h_modes,2);
            num_loadcases = obj.number_of_loadcases;
            dofs = obj.degrees_of_freedom;

            perturbation_displacements = zeros(dofs,num_h_modes,num_loadcases);
            file_path = obj.data_dir + "\" + obj.perturbation_name + "_";
            create_parallel_pool("current");
            parfor iMode = 1:num_h_modes
                file_name =  file_path + h_modes(iMode);
                perturbation_data = load(file_name,"perturbation_disp");
                perturbation_displacements(:,iMode,:) = perturbation_data.perturbation_disp;
            end

           
        end
        %-----
        function obj = combine_data(obj,obj_two)
            if obj.degrees_of_freedom ~= obj_two.degrees_of_freedom
                error("Incompatable perturbations")
            end
 
            if obj.data_dir == obj_two.data_dir
                error("Cannot add a perturbation point to itself")
            end
            
            destination_path = obj.data_dir;
            current_path = obj_two.data_dir;
            num_old_loadcases = obj.number_of_loadcases;
            num_added_loadcases = obj_two.number_of_loadcases;
            
            % perturbation_path = current_path + "\" + obj_two.perturbation_name;
            % file_type_two = obj_two.file_type;
            % 
            % new_perturbation_path = destination_path + "\" + obj.perturbation_name;
            % file_type_one = obj.file_type;

            num_modes = obj_two.num_validation_modes;

            for iMode = 1:num_modes
                perturbation_disp_one = obj.get_displacement(iMode);
                perturbation_disp_two = obj_two.get_displacement(iMode);
                perturbation_disp = cat(3,perturbation_disp_one,perturbation_disp_two);
                obj.save_data(perturbation_disp,iMode)
            end

             obj.num_validation_modes = size(unique([obj.perturbed_modes,obj_two.perturbed_modes]),2);
           
    

            obj.number_of_loadcases = num_old_loadcases + num_added_loadcases;
        end
        %------
        function save_data(obj,perturbation_displacements,validation_modes)
            if nargin == 2
                validation_modes = obj.perturbed_modes;
            end
            num_modes = size(validation_modes,2);
            if isempty(perturbation_displacements)
                return
            end
            file_path =  obj.data_dir + "\" + obj.perturbation_name + "_";
            create_parallel_pool("current");
            parfor iMode = 1:num_modes
                perturbation_disp = squeeze(perturbation_displacements(:,iMode,:));
                file_name = file_path + validation_modes(iMode);
                save(file_name,"-fromstruct",struct("perturbation_disp",perturbation_disp))
            end 
        end
        %------
        function save(Perturbation_Pointer)
            file_dir = Perturbation_Pointer.data_dir + "\" + Perturbation_Pointer.object_name+"_data";
            save(file_dir,Perturbation_Pointer.object_name);
        end
        %------


        %--------------------------
        % overloading
        %------
        function obj = cat(dim,varargin)
            if dim ~= 3
                error("Not implemented")
            end

            num_arrays = size(varargin,2);
            obj = [];
            
            for iArray = 1:num_arrays
                if class(obj) ~= "Perturbation_Pointer"
                    obj = varargin{iArray};
                    continue
                end
                
                next_obj = varargin{iArray};
                if class(next_obj) == "Perturbation_Pointer"
                    obj = obj.combine_data(next_obj);
                end
            end
            obj.save;
        end
        %------
        function obj_size = size(obj,dim)
            obj_size = [obj.degrees_of_freedom,obj.num_validation_modes,obj.number_of_loadcases];
            if nargin == 2
               obj_size = obj_size(dim); 
            end
        end
        %------

    end

end