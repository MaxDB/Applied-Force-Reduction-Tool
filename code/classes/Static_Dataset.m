classdef Static_Dataset
    %Stores the results of applying different loadcases to the system
    properties
        Model
        
        reduced_displacement
        physical_displacement
        restoring_force
        potential_energy
        
        additional_data_type
        tangent_stiffness
        perturbation_displacement
        perturbation_scale_factor

        low_frequency_stiffness
        low_frequency_coupling_gradient
        Dynamic_Validation_Data

        scaffold_points
        static_equilibrium_path_id
        unit_sep_ratios

        Verification_Options
        verified_degree
    end
    methods
        function obj = Static_Dataset(Model,Verification_Opts,varargin)
            

            obj.additional_data_type = Model.Static_Options.additional_data;

            if obj.additional_data_type == "perturbation"
                if isstring(Model.Static_Options.perturbation_scale_factor) && Model.Static_Options.perturbation_scale_factor == "auto"
                    perturbation_sf = select_perturbation_scale_factor(Model);
                else
                    perturbation_sf = Model.Static_Options.perturbation_scale_factor;
                    if isscalar(perturbation_sf)
                        num_h_modes = size(Model.reduced_modes,2) + size(Model.low_frequency_modes,2);
                        perturbation_sf = repmat(perturbation_sf,1,num_h_modes);
                    end
                end
                obj.perturbation_scale_factor = perturbation_sf;
                Model.Static_Options.perturbation_scale_factor = perturbation_sf;
            end
            
            obj.Model = Model;
            data_path = get_data_path(obj);
            if isfolder(data_path)
                rmdir(data_path,"s")
            end



            obj = update_verification_opts(obj,Verification_Opts);

            if nargin == 2
                obj = obj.create_dataset;
                static_dataset_verificiation_plot(obj)
            else
                Properties_Data = varargin{1,1};
                properties = fields(Properties_Data);
                num_properties = size(properties,1);
                for iProp = 1:num_properties
                    property = properties{iProp};
                    obj.(property) = Properties_Data.(property);
                end

            end
        end
        %-----------------------------------------------------------------%
        function obj = create_dataset(obj)
            rom_data_time_start = tic;
            
            
            obj = obj.create_scaffold;
            obj = obj.verify_dataset;
            obj = obj.add_perturbation_data;

            rom_data_time = toc(rom_data_time_start);
            log_message = sprintf("ROM Dataset Created: %.1f seconds" ,rom_data_time);
            logger(log_message,1)
        end
        %-----------------------------------------------------------------%
        function obj = verify_dataset(obj)
            %add loadcases until convergence
            verification_algorithm = obj.Verification_Options.verification_algorithm;
            verification_time_start = tic;
            switch verification_algorithm
                case "sep_to_edge"
                    obj = sep_verification(obj);
                case "sep_from_origin"
                    obj = sep_from_origin_verification(obj);
                case "sep_grow"
                    obj = sep_grow_verification(obj);
            end
            % switch additional_data_type
            %     case "stiffness"
            %         obj = minimum_stiffness_degree(obj);
            % end
            verification_time = toc(verification_time_start);
            log_message = sprintf("ROM Verified: %.1f seconds" ,verification_time);
            logger(log_message,1)
        end
        %-----------------------------------------------------------------%
        function obj = add_additional_data(obj,Static_Opts)
            if obj.Model.Static_Options.additional_data ~= "none"
                error("adding to existing data not implemented")
            end
            obj.Model = obj.Model.update_static_opts(Static_Opts);

            obj.additional_data_type = Static_Opts.additional_data;
            if obj.additional_data_type == "perturbation"
                if isstring(obj.Model.Static_Options.perturbation_scale_factor) && obj.Model.Static_Options.perturbation_scale_factor == "auto"
                    perturbation_sf = select_perturbation_scale_factor(obj.Model);
                else
                    perturbation_sf = obj.Model.Static_Options.perturbation_scale_factor;
                    if isscalar(perturbation_sf)
                        num_h_modes = size(obj.Model.reduced_modes,2) + size(obj.Model.low_frequency_modes,2);
                        perturbation_sf = repmat(perturbation_sf,1,num_h_modes);
                    end
                end
                obj.perturbation_scale_factor = perturbation_sf;
                obj.Model.Static_Options.perturbation_scale_factor = perturbation_sf;
            end
            
            Closest_Point.initial_disp = obj.get_dataset_values("physical_displacement");
            Closest_Point.initial_force = obj.get_dataset_values("restoring_force");
            
            dummy_loads = nan(size(Closest_Point.initial_force));

            [~,~,~,~,additional_data] = obj.Model.add_point(dummy_loads,obj.additional_data_type,Closest_Point);
            obj = obj.update_additional_data(additional_data);
        end
        %-----------------------------------------------------------------%
        function obj = add_perturbation_data(obj)
            if obj.additional_data_type ~= "stiffness"
                return
            end
            perturbation_time_start = tic;
            Static_Opts = obj.Model.Static_Options;

            if nargin == 2
                if isempty(L_modes) || L_modes == 0
                    return
                end
            elseif nargin == 1
                L_modes = 1:Static_Opts.num_validation_modes;
            end


            r_modes = obj.Model.reduced_modes;
            L_modes(ismember(L_modes,r_modes)) = [];

            found_L_modes = obj.Model.low_frequency_modes;
            unknown_L_modes = setdiff(L_modes,found_L_modes);
            calc_eigenproblem = ~isempty(unknown_L_modes);

            if calc_eigenproblem
                mass = obj.Model.mass;
                stiffness = obj.Model.stiffness;

                [full_evecs,full_evals] = eigs(stiffness,mass,max(L_modes),"smallestabs");
                full_evals = diag(full_evals);

                new_L_modes = 1:max(L_modes);
                new_L_modes(ismember(new_L_modes,r_modes)) = [];

                obj.Model.low_frequency_modes = new_L_modes;
                obj.Model.low_frequency_eigenvalues = full_evals(new_L_modes);

                evec_L =  full_evecs(:,new_L_modes);
                matrix_data = whos("evec_L");
                if matrix_data.bytes/1024 > obj.Model.Static_Options.max_matrix_size
                    data_path = obj.Model.get_data_path;
                    evec_L = Large_Matrix_Pointer(evec_L,data_path,"validation_eigenvectors");
                end
                obj.Model.low_frequency_eigenvectors = evec_L;
            end

            obj = apply_small_force(obj,L_modes);
            perturbation_time = toc(perturbation_time_start);
            log_message = sprintf("Perturbation data created: %.1f seconds" ,perturbation_time);
            logger(log_message,2)
        end
        %-----------------------------------------------------------------%
        function obj = add_validation_data(obj,validation_modes)


            if isempty(obj.get_dataset_values("perturbation_displacement"))
                error("Cannot find perturbation data")
            end


            r_modes = obj.Model.reduced_modes;
            validation_modes(ismember(validation_modes,r_modes)) = [];

            found_L_modes = obj.Model.low_frequency_modes;
            unknown_L_modes = setdiff(validation_modes,found_L_modes);
            if ~isempty(unknown_L_modes)
                switch obj.additional_data_type
                    case "perturbation"
                        error("Perturbations for mode(s) " + join(string(unknown_L_modes),", ") + " are not in static dataset")
                    case "stiffness"
                        error("not implemented")
                end
            end
            if ~isempty(obj.Dynamic_Validation_Data)
                if isequal(obj.Dynamic_Validation_Data.current_L_modes,validation_modes)
                    return
                end
            end

            validation_data_start = tic;
            [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_perturbation_data(obj,validation_modes);

            validation_data_time = toc(validation_data_start);
            log_message = sprintf("Processed validation data: %.1f seconds" ,validation_data_time);
            logger(log_message,3)

            obj.low_frequency_stiffness = h_stiffness;
            obj.low_frequency_coupling_gradient = h_coupling_gradient;

            obj.Dynamic_Validation_Data.current_L_modes = validation_modes;
            obj.Dynamic_Validation_Data.h_stiffness_0 = h_stiffness_0;
            obj.Dynamic_Validation_Data.h_coupling_gradient_0 = h_coupling_gradient_0;

            obj.verified_degree = repmat(obj.verified_degree,1,2);
            % minimum_degree_start = tic;
            % % obj = minimum_h_degree(obj);
            %
            %
            % minimum_degree_time = toc(minimum_degree_start);
            % log_message = sprintf("Validating validation polynomials: %.1f seconds" ,minimum_degree_time);
            % logger(log_message,3)

            validation_dataset_verificiation_plot(obj)

        end
        %-----------------------------------------------------------------%
        function obj = create_scaffold(obj,varargin)
            %-------------------------------------------------------------------------%
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            sep_density = 1;
            limit_scale_factor = 1;
            verification_algorithm = obj.Verification_Options.verification_algorithm;

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "density"
                        sep_density = keyword_values{arg_counter};
                    case "sf"
                        limit_scale_factor = keyword_values{arg_counter};
                    case "algorithm"
                        verification_algorithm = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------%
            rom_scaffold_time_start = tic;
            
            switch verification_algorithm
                case {"sep_to_edge","sep_from_origin"}
                    r_modes = obj.Model.reduced_modes;
                    num_modes = length(r_modes);

                    full_unit_force_ratios = add_sep_ratios(num_modes,sep_density);

                    found_sep_ratios = obj.unit_sep_ratios;
                    unit_force_ratios = add_sep_ratios(num_modes,sep_density,found_sep_ratios);

                    member_map = ismembertol(full_unit_force_ratios',unit_force_ratios',"ByRows",true)';
                    sep_map = [find(~member_map),find(member_map)];

                    scaled_force_ratios = scale_sep_ratios(unit_force_ratios,obj.Model.calibrated_forces);
                    scaled_force_ratios = scaled_force_ratios*limit_scale_factor;

                    [r,theta,f,E,sep_id,additional_data] = obj.Model.add_sep(scaled_force_ratios,obj.additional_data_type,1);

                    obj = obj.update_data(r,theta,f,E,sep_id,additional_data,unit_force_ratios);
                    scaf_points = obj.scaffold_points == 1;
                    obj.static_equilibrium_path_id(scaf_points) = sep_map(obj.static_equilibrium_path_id(scaf_points));
                    obj.unit_sep_ratios(:,sep_map) = obj.unit_sep_ratios;
                case "sep_grow"
                    if num_args == 0
                        return
                    end

                    % obj = obj.create_scaffold("algorithm","sep_from_origin",varargin{:});
                    r_modes = obj.Model.reduced_modes;
                    num_modes = length(r_modes);

                    full_unit_force_ratios = add_sep_ratios(num_modes,sep_density);

                    unit_force_ratios = add_sep_ratios(num_modes,sep_density,[]);

                    member_map = ismembertol(full_unit_force_ratios',unit_force_ratios',"ByRows",true)';
                    sep_map = [find(~member_map),find(member_map)];
                    

                    scaled_force_ratios = scale_sep_ratios(unit_force_ratios,obj.Model.calibrated_forces);
                    scaled_force_ratios = scaled_force_ratios*limit_scale_factor;
                    
                    [r,theta,f,E,additional_data] = obj.Model.add_point(scaled_force_ratios,obj.additional_data_type);
                    sep_id = 1:size(r,2);

                    obj = obj.update_data(r,theta,f,E,sep_id,additional_data,unit_force_ratios);
                    scaf_points = obj.scaffold_points == 1;
                    num_scaf_points = nnz(scaf_points);
                    sep_map_rep = repmat(sep_map,num_scaf_points/size(sep_map,2));

                    % obj.static_equilibrium_path_id(scaf_points) = sep_map_rep(obj.static_equilibrium_path_id(scaf_points));
                    % obj.unit_sep_ratios = obj.unit_sep_ratios(:,sep_map);

            end
            rom_scaffold_time = toc(rom_scaffold_time_start);
            log_message = sprintf("Dataset scaffold created: %.1f seconds" ,rom_scaffold_time);
            logger(log_message,2)


        end
        %-----------------------------------------------------------------%
        function obj = update_data(obj,r,x,f,V,sep_id,additional_data,found_force_ratios)
            obj.reduced_displacement = [obj.get_dataset_values("reduced_displacement"),r];
            obj.physical_displacement = [obj.get_dataset_values("physical_displacement"),x];
            obj.restoring_force = [obj.get_dataset_values("restoring_force"),f];
            obj.potential_energy = [obj.get_dataset_values("potential_energy"),V];
            
            sep_id_shift = size(obj.unit_sep_ratios,2);
            % sep_id_shift = max(obj.static_equilibrium_path_id);
            obj.static_equilibrium_path_id = [obj.static_equilibrium_path_id,sep_id+sep_id_shift];

            num_points = size(V,2);
            if nargin == 8
                obj.unit_sep_ratios = [obj.unit_sep_ratios,found_force_ratios];
                obj.scaffold_points = [obj.scaffold_points,true(1,num_points)];
            else
                obj.scaffold_points = [obj.scaffold_points,false(1,num_points)];
            end
            
            obj = update_additional_data(obj,additional_data);
        end
        %-----------------------------------------------------------------%
        function obj = update_additional_data(obj,additional_data)
            switch obj.additional_data_type
                case "stiffness"
                    current_stiffness = obj.get_dataset_values("tangent_stiffness");
                    if isempty(current_stiffness)
                       data_path = get_data_path(obj);
                       current_stiffness = Sparse_Stiffness_Pointer(obj.Model,"path",data_path + "stiffness");
                    end
                    obj.tangent_stiffness = cat(3,current_stiffness,additional_data);
                    
                case "perturbation"
                    obj.perturbation_displacement = cat(3,obj.get_dataset_values("perturbation_displacement"),additional_data);
            end
        end
        %-----------------------------------------------------------------%
        function obj = update_model(obj,added_modes,Static_Opts,Calibration_Opts)
            TEMP_PROPERTY_VALUE = "unloaded";

            system_name = obj.Model.system_name;
            energy_limit = obj.Model.energy_limit;
            initial_modes = sort([obj.Model.reduced_modes,added_modes],"ascend");
            old_mode_map = ismember(initial_modes,obj.Model.reduced_modes);
            old_L_modes = obj.Model.low_frequency_modes;
            old_h_modes = [obj.Model.reduced_modes,old_L_modes];

            if nargin > 2
                Static_Opts.static_solver = obj.Model.Static_Options.static_solver;
                Static_Opts.additional_data = obj.Model.Static_Options.additional_data;
                Static_Opts.num_validation_modes = obj.Model.Static_Options.num_validation_modes;
                obj.Model = obj.Model.update_static_opts(Static_Opts);
            end
            Static_Opts = obj.Model.Static_Options;

            if nargin < 4
                Calibration_Opts = obj.Model.Calibration_Options;
            end
            %-------------------------------------------------%
            %load old data
            static_data_properties = properties(obj);
            num_properties = size(static_data_properties,1);
            for iProp = 1:num_properties
                static_data_property = static_data_properties{iProp,1};
                if isstring(obj.(static_data_property)) && obj.(static_data_property) == TEMP_PROPERTY_VALUE
                    obj.(static_data_property) = obj.get_dataset_values(static_data_property);
                end
            end
            %-------------------------------------------------%

            obj.Model = Dynamic_System(system_name,energy_limit,initial_modes,"calibration_opts",Calibration_Opts,"static_opts",Static_Opts);


            %-------------------------------------------------%
            num_modes = size(initial_modes,2);
            

            old_sep_ratios = obj.unit_sep_ratios;
            num_seps = size(old_sep_ratios,2);
            new_sep_ratios = zeros(num_modes,num_seps);
            new_sep_ratios(old_mode_map,:) = old_sep_ratios;
            obj.unit_sep_ratios = new_sep_ratios;

            r_evecs = obj.Model.reduced_eigenvectors;
            mass = obj.Model.mass;
            obj.reduced_displacement = r_evecs'*mass*obj.get_dataset_values("physical_displacement");

            old_force = obj.get_dataset_values("restoring_force");
            num_loadcases = size(old_force,2);
            new_force = zeros(num_modes,num_loadcases);
            new_force(old_mode_map,:) = old_force;
            obj.restoring_force = new_force;

            if obj.additional_data_type == "perturbation"
                L_modes = obj.Model.low_frequency_modes;
                h_modes = [obj.Model.reduced_modes,L_modes];
                h_map = arrayfun(@(h_mode) find(old_h_modes == h_mode),h_modes);
                perturbation_disp = obj.get_dataset_values("perturbation_displacement");
                obj.perturbation_displacement = perturbation_disp(:,h_map,:);
            end
            %-------------------------------------------------%
           
        end
        %-----------------------------------------------------------------%


        %-----------------------------------------------------------------%
        function obj = update_verification_opts(obj,Verification_Opts)
            % Default static options
            Default_Verification_Opts = read_default_options("verification");         
            obj.Verification_Options = update_options(Default_Verification_Opts,obj.Verification_Options,Verification_Opts);
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %Overloading
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function sz = size(obj,dim)
            sz = size(obj.restoring_force);
            if nargin == 2
                sz = sz(dim);
            end
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %Helpers
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function save_data(Static_Data)
            SEPERATELY_SAVED_PROPERTIES = [
                "physical_displacement", ...
                "low_frequency_stiffness", ...
                "low_frequency_coupling_gradient"];

            data_path = get_data_path(Static_Data);
            system_data_path = split(data_path,"\");
            system_data_path = join(system_data_path(1:(end-2)),"\") + "\";
            % if exist(system_data_path,"dir")
            %     % check if any data is unloaded
            %     for iProp = 1:size(SEPERATELY_SAVED_PROPERTIES,2)
            %         saved_property = SEPERATELY_SAVED_PROPERTIES(iProp);
            %         if isstring(Static_Data.(saved_property)) && Static_Data.(saved_property) == "unloaded"
            %             Static_Data.(saved_property) = Static_Data.get_dataset_values(saved_property);
            %         end
            %     end
            %     rmdir(system_data_path,"s")
            % end
            % mkdir(data_path)

            if ~isfolder(system_data_path)
                mkdir(data_path)
            end
            
            num_properties = size(SEPERATELY_SAVED_PROPERTIES,2);
            for iProperty = 1:num_properties
                Static_Data = save_property_data(SEPERATELY_SAVED_PROPERTIES(iProperty),Static_Data);
            end

            
            save(data_path + "Static_Data.mat","Static_Data","-v7.3")

            function Static_Data = save_property_data(static_data_property,Static_Data)
                TEMP_PROPERTY_VALUE = "unloaded";
                if isempty(Static_Data.(static_data_property))
                    return
                end
                value = Static_Data.(static_data_property);
                if isstring(value) && value == "unloaded"
                    return %dont overwrite data with "unloadaed"
                end
                Static_Data.(static_data_property) = TEMP_PROPERTY_VALUE;
                
                file_path = get_data_path(Static_Data) + static_data_property + ".mat";
                save(file_path,"value","-v7.3")
            end
        end
        %-----------------------------------------------------------------%
        function value = get_dataset_values(obj,static_data_property)
            if ~isstring(obj.(static_data_property))
                value = obj.(static_data_property);
                return
            end
            file_path = get_data_path(obj) + static_data_property + ".mat";
            if ~isfile(file_path)
                error("Can't locate " + static_data_property)
            end
            static_property_data = load(file_path);
            value = static_property_data.value;
        end
        %-----------------------------------------------------------------%
        function data_path = get_data_path(obj)
            r_modes = obj.Model.reduced_modes;
            mode_id = join(string(r_modes),"");
            data_path = "data\" + obj.Model.system_name + "_" + mode_id + "\static_data\";
        end
        %-----------------------------------------------------------------%
        function[Static_Data,Static_Data_Removed] = remove_loadcases(Static_Data,removal_index)
            Static_Data_Removed = Static_Data;

            Static_Data_Removed.reduced_displacement = Static_Data_Removed.reduced_displacement(:,removal_index);
            Static_Data_Removed.physical_displacement = Static_Data_Removed.physical_displacement(:,removal_index);
            Static_Data_Removed.potential_energy = Static_Data_Removed.potential_energy(:,removal_index);
            Static_Data_Removed.restoring_force = Static_Data_Removed.restoring_force(:,removal_index);
            
            Static_Data_Removed.scaffold_points = Static_Data_Removed.scaffold_points(:,removal_index);
            Static_Data_Removed.static_equilibrium_path_id = Static_Data_Removed.static_equilibrium_path_id(:,removal_index);

            switch Static_Data_Removed.additional_data_type
                case "stiffness"
                    Static_Data_Removed.tangent_stiffness = Static_Data_Removed.tangent_stiffness(:,:,removal_index);
                case "perturbation"
                    Static_Data_Removed.perturbation_displacement = Static_Data_Removed.perturbation_displacement(:,:,removal_index);
            end
            %---------------------------------------------------------------------------------------%
            Static_Data.reduced_displacement = Static_Data.reduced_displacement(:,~removal_index);
            Static_Data.physical_displacement = Static_Data.physical_displacement(:,~removal_index);
            Static_Data.potential_energy = Static_Data.potential_energy(:,~removal_index);
            Static_Data.restoring_force = Static_Data.restoring_force(:,~removal_index);
            
            Static_Data.scaffold_points = Static_Data.scaffold_points(:,~removal_index);
            Static_Data.static_equilibrium_path_id = Static_Data.static_equilibrium_path_id(:,~removal_index);

            switch Static_Data.additional_data_type
                case "stiffness"
                    Static_Data.tangent_stiffness = Static_Data.tangent_stiffness(:,:,~removal_index);
                case "perturbation"
                    Static_Data.perturbation_displacement = Static_Data.perturbation_displacement(:,:,~removal_index);
            end
        end
        %-----------------------------------------------------------------%
        function [loadcases_found,r,x,f,E,additional_data] = contains_loadcase(obj,loadcases)
            MAX_DIFF = 0.05;
            
            found_force = obj.get_dataset_values("restoring_force");
            num_r_modes = size(found_force,1);
            num_dofs = obj.Model.num_dof;
            num_checked_loadcases = size(loadcases,2);
            loadcases_found = false(1,num_checked_loadcases);
            
            matched_counter = 0;
            r = zeros(num_r_modes,0);
            x = zeros(num_dofs,0);
            f = zeros(num_r_modes,0);
            E = zeros(1,0);
            

            r_old = obj.get_dataset_values("reduced_displacement");
            x_old = obj.get_dataset_values("physical_displacement");
            f_old = obj.get_dataset_values("restoring_force");
            E_old = obj.get_dataset_values("potential_energy");

            switch obj.additional_data_type
                case "stiffness"
                    tangent_stiff_old = obj.get_dataset_values("tangent_stiffness");
                    if issparse(tangent_stiff_old)
                        additional_data = sparse(num_dofs,num_dofs,0);
                    else
                        additional_data = zeros(num_dofs,num_dofs,0);
                    end
                case "perturbation"
                    perturbation_disp_old = obj.get_dataset_values("perturbation_displacement");
                    num_h_modes = size(perturbation_disp_old,2);
                    additional_data = zeros(num_dofs,num_h_modes,0);
                case "none"
                    additional_data = [];
            end

            for iLoad = 1:num_checked_loadcases
                loadcase = loadcases(:,iLoad);
                bounds = abs(found_force./loadcase - 1);
                bounds(isnan(bounds)) = 0;
                matching_index = all(bounds < MAX_DIFF,1);
                if any(matching_index)
                    matched_counter = matched_counter + 1;
                    loadcases_found(iLoad) = true;
                    loadcase_index = find(matching_index); %problematic if two are detected, should pick closest (but also bad if happens)
                    if size(loadcase_index,2) > 1
                        [~,min_index] = min(sqrt(sum((found_force(:,loadcase_index)-loadcase).^2,1)));
                        loadcase_index = loadcase_index(min_index);
                    end
                    r(:,matched_counter) = r_old(:,loadcase_index);
                    x(:,matched_counter) = x_old(:,loadcase_index);
                    f(:,matched_counter) = f_old(:,loadcase_index);
                    E(:,matched_counter) = E_old(:,loadcase_index);

                    switch obj.additional_data_type
                        case "stiffness"
                            additional_data(:,:,matched_counter) = tangent_stiff_old(:,:,loadcase_index);
                        case "perturbation"
                            additional_data(:,:,matched_counter) = perturbation_disp_old(:,:,loadcase_index);
                        case "none"
                            additional_data = [];
                    end

                end
            end

        end
        %-----------------------------------------------------------------%
        function [h_modes,h_eval,h_evec] = get_current_h_data(obj)
            r_modes = obj.Model.reduced_modes;
            r_eval = obj.Model.reduced_eigenvalues;
            r_evec = obj.Model.reduced_eigenvectors;

            current_L_modes = obj.Dynamic_Validation_Data.current_L_modes;
            L_modes = obj.Model.low_frequency_modes;
            L_map = ismember(L_modes,current_L_modes);
            L_eval = obj.Model.low_frequency_eigenvalues(L_map);
            L_evec = obj.Model.low_frequency_eigenvectors(:,L_map);
            
            h_modes = [r_modes,current_L_modes];
            h_eval = [r_eval;L_eval];
            h_evec = [r_evec,L_evec];
        end
        %-----------------------------------------------------------------%
    end
end