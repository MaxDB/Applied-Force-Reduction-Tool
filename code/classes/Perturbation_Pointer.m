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

            obj.number_of_loadcases = size(perturbation_disp,3);
            obj.num_validation_modes = size(perturbation_disp,2);
            obj.degrees_of_freedom = size(perturbation_disp,1);

            obj.perturbed_modes = [Model.reduced_modes,Model.low_frequency_modes];


            obj.save_data(perturbation_disp);
            obj.save();
        end
        %-----
        function perturbation_displacements = get_displacement(obj,h_modes)
            num_h_modes = size(h_modes,2);
            num_loadcases = obj.number_of_loadcases;
            dofs = obj.degrees_of_freedom;

            perturbation_displacements = zeros(dofs,num_h_modes,num_loadcases);
            file_path = obj.data_dir + "\" + obj.perturbation_name + "_";
            parfor iMode = 1:num_h_modes
                file_name =  file_path + h_modes(iMode);
                perturbation_data = load(file_name,"perturbation_disp");
                perturbation_displacements(:,iMode,:) = perturbation_data.perturbation_disp;
            end

           
        end

        %-----
        function save_data(obj,perturbation_displacements)
            validation_modes = obj.perturbed_modes;
            num_modes = size(validation_modes,2);

            file_path =  obj.data_dir + "\" + obj.perturbation_name + "_";
            % parfor iMode = 1:num_modes
            for iMode = 1:num_modes
                perturbation_disp = squeeze(perturbation_displacements(:,iMode,:));
                file_name = file_path + validation_modes(iMode);
                save(file_name,"perturbation_disp")
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
        function obj_size = size(obj,dim)
            obj_size = [obj.degrees_of_freedom,obj.num_validation_modes,obj.number_of_loadcases];
            if nargin == 2
               obj_size = obj_size(dim); 
            end
        end
        %------

    end

end