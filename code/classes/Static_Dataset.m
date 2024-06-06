classdef Static_Dataset
    %Stores the results of applying different loadcases to the system
    properties
        Model
        
        reduced_displacement
        condensed_displacement
        restoring_force
        potential_energy
        
        additional_data_type
        tangent_stiffness
        perturbation_displacement
        perturbation_scale_factor

        low_frequency_stiffness
        low_frequency_coupling_gradient
        Dynamic_Validation_Data

        static_equilibrium_path_id
        unit_sep_ratios

        Validation_Options
        validated_degree
    end
    methods
        function obj = Static_Dataset(Model,Validation_Opts)
            rom_data_time_start = tic;
            obj.Model = Model;
            obj.additional_data_type = Model.Static_Options.additional_data;
            obj.perturbation_scale_factor = Model.Static_Options.perturbation_scale_factor;

            obj = update_validation_opts(obj,Validation_Opts);

            obj = obj.create_dataset;
        
            
            rom_data_time = toc(rom_data_time_start);
            log_message = sprintf("ROM Dataset Created: %.1f seconds" ,rom_data_time);
            logger(log_message,1)
        end
        %-----------------------------------------------------------------%
        function obj = create_dataset(obj)
            %using hyper sphere but may update
            rom_scaffold_time_start = tic;

            obj = obj.create_scaffold(obj.additional_data_type);

            rom_scaffold_time = toc(rom_scaffold_time_start);
            log_message = sprintf("Dataset scaffold created: %.1f seconds" ,rom_scaffold_time);
            logger(log_message,2)

            %add loadcases until convergence
            validation_time_start = tic;
            switch obj.Validation_Options.validation_algorithm
                case "sep"
                    obj = sep_validation(obj);
                case "sep_points"
                    obj = sep_points_validation(obj);
                case "grid"
                    obj = grid_validation(obj);
            end
            % switch additional_data_type
            %     case "stiffness"
            %         obj = minimum_stiffness_degree(obj);
            % end
            validation_time = toc(validation_time_start);
            log_message = sprintf("ROM Dataset Validated: %.1f seconds" ,validation_time);
            logger(log_message,1)
        end
        %-----------------------------------------------------------------%
        function obj = add_validation_data(obj,L_modes)
            switch obj.additional_data_type
                case "stiffness"
                    r_modes = obj.Model.reduced_modes;
                    L_modes(ismember(L_modes,r_modes)) = [];

                    Validation_Data = obj.Dynamic_Validation_Data;
                    calc_eigenproblem = isempty(Validation_Data);
                    if ~calc_eigenproblem
                        if isequal(Validation_Data.current_L_modes,L_modes)
                            return
                        end

                        if Validation_Data.largest_L_mode < max(L_modes)
                            calc_eigenproblem = 1;
                        end
                    end
                    if calc_eigenproblem
                        mass = obj.Model.mass;
                        stiffness = obj.Model.stiffness;

                        [full_evecs,full_evals] = eigs(stiffness,mass,max(L_modes),"smallestabs");
                        full_evals = diag(full_evals);


                        obj.Dynamic_Validation_Data.full_eigenvectors = full_evecs;
                        obj.Dynamic_Validation_Data.full_eigenvalues = full_evals;
                        obj.Dynamic_Validation_Data.largest_L_mode = max(L_modes);
                    end
                    obj.Dynamic_Validation_Data.current_L_modes = L_modes;
                    L_evecs = obj.Dynamic_Validation_Data.full_eigenvectors(:,L_modes);

                    [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_h_error_data(obj,L_evecs);

                    obj.low_frequency_stiffness = h_stiffness;
                    obj.low_frequency_coupling_gradient = h_coupling_gradient;

                    obj.Dynamic_Validation_Data.h_stiffness_0 = h_stiffness_0;
                    obj.Dynamic_Validation_Data.h_coupling_gradient_0 = h_coupling_gradient_0;

                    obj = minimum_h_degree(obj);
                case "perturbation"
                    r_modes = obj.Model.reduced_modes;
                    L_modes(ismember(L_modes,r_modes)) = [];

                    found_L_modes = obj.Model.low_frequency_modes;
                    unknown_L_modes = setdiff(L_modes,found_L_modes);
                    if ~isempty(unknown_L_modes)
                        error("Perturbations for mode(s) " + join(string(unknown_L_modes),", ") + " are not in static dataset")
                    end
                    if ~isempty(obj.Dynamic_Validation_Data)
                        if isequal(obj.Dynamic_Validation_Data.current_L_modes,L_modes)
                            return
                        end
                    end


                    [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_perturbation_data(obj,L_modes);

                    obj.low_frequency_stiffness = h_stiffness;
                    obj.low_frequency_coupling_gradient = h_coupling_gradient;

                    obj.Dynamic_Validation_Data.current_L_modes = L_modes;
                    obj.Dynamic_Validation_Data.h_stiffness_0 = h_stiffness_0;
                    obj.Dynamic_Validation_Data.h_coupling_gradient_0 = h_coupling_gradient_0;

                    obj = minimum_h_degree(obj);

            end
        end
        %-----------------------------------------------------------------%
        function obj = create_scaffold(obj,additional_data_type)

            r_modes = obj.Model.reduced_modes;
            num_modes = length(r_modes);
            
            found_sep_ratios = obj.unit_sep_ratios;
            unit_force_ratios = add_sep_ratios(num_modes,1,found_sep_ratios);
            scaled_force_ratios = scale_sep_ratios(unit_force_ratios,obj.Model.calibrated_forces);

      
            
            
            [r,theta,f,E,sep_id,additional_data] = obj.Model.add_sep(scaled_force_ratios,additional_data_type,1);
            
            obj = obj.update_data(r,theta,f,E,sep_id,additional_data,unit_force_ratios);

            % obj.unit_sep_ratios = [found_sep_ratios,unit_force_ratios];
            % obj.reduced_displacement = r;
            % obj.condensed_displacement = theta;
            % obj.restoring_force = f;
            % obj.potential_energy = E;
            % obj.static_equilibrium_path_id = sep_id;
            
            % switch additional_data_type
            %     case "stiffness"
            %         obj.tangent_stiffness = additional_data;
            %     case "perturbation"
            %         obj.perturbation_displacement = additional_data;
            % end
            

        end
        %-----------------------------------------------------------------%
        function obj = update_data(obj,r,theta,f,V,sep_id,additional_data,found_force_ratios)
            obj.reduced_displacement = [obj.reduced_displacement,r];
            obj.condensed_displacement = [obj.condensed_displacement,theta];
            obj.restoring_force = [obj.restoring_force,f];
            obj.potential_energy = [obj.potential_energy,V];
            
            sep_id_shift = size(obj.unit_sep_ratios,2);
            % sep_id_shift = max(obj.static_equilibrium_path_id);
            obj.static_equilibrium_path_id = [obj.static_equilibrium_path_id,sep_id+sep_id_shift];

            if nargin == 8
                obj.unit_sep_ratios = [obj.unit_sep_ratios,found_force_ratios];
            end

            switch obj.additional_data_type
                case "stiffness"
                    obj.tangent_stiffness = cat(3,obj.tangent_stiffness,additional_data);
                case "perturbation"
                    obj.perturbation_displacement = cat(3,obj.perturbation_displacement,additional_data);
            end
        end
        %-----------------------------------------------------------------%
        function obj = update_model(obj,added_modes,Calibration_Opts,Static_Opts)
            system_name = obj.Model.system_name;
            energy_limit = obj.Model.energy_limit;
            initial_modes = sort([obj.Model.reduced_modes,added_modes],"ascend");
            old_mode_map = ismember(initial_modes,obj.Model.reduced_modes);
           
            Static_Opts.static_solver = obj.Model.Static_Options.static_solver;
            Static_Opts.additional_data = obj.Model.Static_Options.additional_data;
            Static_Opts.num_validation_modes = obj.Model.Static_Options.num_validation_modes;

            obj.Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);


            %-------------------------------------------------%
            num_modes = size(initial_modes,2);
            

            old_sep_ratios = obj.unit_sep_ratios;
            num_seps = size(old_sep_ratios,2);
            new_sep_ratios = zeros(num_modes,num_seps);
            new_sep_ratios(old_mode_map,:) = old_sep_ratios;
            obj.unit_sep_ratios = new_sep_ratios;

            r_evecs = obj.Model.reduced_eigenvectors;
            mass = obj.Model.mass;
            obj.reduced_displacement = r_evecs'*mass*obj.condensed_displacement;

            old_force = obj.restoring_force;
            num_loadcases = size(old_force,2);
            new_force = zeros(num_modes,num_loadcases);
            new_force(old_mode_map,:) = old_force;
            obj.restoring_force = new_force;
            %-------------------------------------------------%
           
        end
        %-----------------------------------------------------------------%


        %-----------------------------------------------------------------%
        function obj = update_validation_opts(obj,Validation_Opts)
            % Default static options
            Default_Validation_Opts = read_default_options("validation");         
            obj.Validation_Options = update_options(Default_Validation_Opts,obj.Validation_Options,Validation_Opts);
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %Overloading
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function sz = size(obj)
            sz = length(obj.potential_energy);
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %Helpers
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function strongly_coupled_dofs = rank_coupling(obj,min_rating)
            r = obj.reduced_displacement;
            r_evec = obj.Model.reduced_eigenvectors;
            theta = obj.condensed_displacement - r_evec*r;
            sep_id = obj.static_equilibrium_path_id;
            M = obj.Model.mass;

            num_seps = max(sep_id);
            sep_ends = arrayfun(@(id) find(sep_id == id,1,"last"),1:num_seps);

            theta_sep_ends = theta(:,sep_ends);

            couple_rating = zeros(size(theta_sep_ends));
            for iSep = 1:num_seps
                theta_sep_end = theta_sep_ends(:,iSep);
                mass_theta_prod = M*theta_sep_end;
                disp_norm = theta_sep_end'*mass_theta_prod;
                theta_contribution = theta_sep_end.*mass_theta_prod;

                couple_rating(:,iSep) = abs(theta_contribution/disp_norm);
            end
            
            couple_rating = max(couple_rating,[],2);
            normalised_couple_rating = couple_rating/max(couple_rating);
            strongly_coupled_dofs = find(normalised_couple_rating>min_rating);

        end
        %-----------------------------------------------------------------%
        function save_data(Static_Data)
            data_path = get_data_path(Static_Data);
            if exist(data_path,"dir")
                rmdir(data_path,"s")
            end
            mkdir(data_path)
            save(data_path + "Static_Data.mat","Static_Data","-v7.3")
        end
        %-----------------------------------------------------------------%
        function data_path = get_data_path(obj)
            r_modes = obj.Model.reduced_modes;
            mode_id = join(string(r_modes),"");
            data_path = "data\" + obj.Model.system_name + "_" + mode_id + "\";
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        %Plotting
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function ax = plot_force(obj,ax,r,f,sep_id)
            EQUILIBRIUM = [0;0];
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);

            if ~exist("r","var")
                r = obj.reduced_displacement;
                f = obj.restoring_force;
                sep_id = obj.static_equilibrium_path_id;
            end

            num_seps = max(sep_id);
            
            switch num_modes
                case 1
                    figure
                    ax = axes(gcf);
                    xlabel("r_" + modes(1))
                    ylabel("f_" + modes(1))
                    box on

                    hold on
                    for iSep = 1:num_seps
                        sep_span = sep_id==iSep;
                        r_sep = r(:,sep_span);
                        f_sep = f(:,sep_span);
                        [~,sep_order] = sort(abs(f_sep(1,:)),"ascend");

                        x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                        y = [0,f_sep(1,sep_order)];
                        plot(x,y,'.-')
                    end
                    plot(EQUILIBRIUM(1),0,'k.',"MarkerSize",8)
                    hold off

                case 2
                    figure
                    tiledlayout(1,num_modes)
                    ax = cell(1,num_modes);
                    for iMode = 1:num_modes
                        ax{1,iMode} = nexttile;
                        xlabel("r_" + modes(1))
                        ylabel("r_" + modes(2))
                        zlabel("f_" + modes(iMode))
                        box on

                        hold on
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            r_sep = r(:,sep_span);
                            f_sep = f(:,sep_span);
                            [~,sep_order] = sort(max(abs(f_sep),[],1),"ascend");

                            x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                            y = [EQUILIBRIUM(2),r_sep(2,sep_order)];
                            z = [0,f_sep(iMode,sep_order)];
                            plot3(x,y,z,'.-')
                        end
                        plot3(EQUILIBRIUM(1),EQUILIBRIUM(2),0,'k.',"MarkerSize",8)
                        hold off
                    end
            end
        end
        %-----------------------------------------------------------------%
        function ax = plot_condensed_displacement(obj,plotted_outputs,ax,r,theta,sep_id)
            EQUILIBRIUM = [0;0];

            num_plots = length(plotted_outputs);
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);

            if ~exist("r","var")
                r = obj.reduced_displacement;
                theta = obj.condensed_displacement;
                sep_id = obj.static_equilibrium_path_id;
                f = obj.restoring_force;
            end

            r_evec = obj.Model.reduced_eigenvectors;
            % theta = theta - r_evec*r;

            num_seps = max(sep_id);

            if ~exist("ax","var")
                figure
                 tiledlayout("flow")
                ax = cell(1,num_plots);
                for iPlot = 1:num_plots
                    
                   
                    ax{1,iPlot} = nexttile;
                end
            end

            if ~iscell(ax)
                ax_temp = ax;
                ax = cell(1,length(ax_temp));
                for iAx = 1:length(ax_temp)
                    ax{1,iAx} = ax_temp(iAx);
                end
                
            end
            

            switch num_modes
                case 1
                    
                    for iPlot = 1:num_plots
                        xlabel(ax{1,iPlot},"r_" + modes)
                        ylabel(ax{1,iPlot},"\theta_{" + plotted_outputs(iPlot) + "}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            r_sep = r(:,sep_span);
                            theta_sep = theta(:,sep_span);
                            [~,sep_order] = sort(abs(f(1,sep_span)),"ascend");
                            
                            x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                            y = [0,theta_sep(plotted_outputs(iPlot),sep_order)];
                            plot(ax{1,iPlot},x,y,'.-')
                        end
                        plot(ax{1,iPlot},EQUILIBRIUM(1),0,'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end
                case 2
             

                    for iPlot = 1:num_plots
                        xlabel(ax{1,iPlot},"r_" + modes(1))
                        ylabel(ax{1,iPlot},"r_" + modes(2))
                        zlabel(ax{1,iPlot},"\theta_{" + plotted_outputs(iPlot) + "}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            r_sep = r(:,sep_span);
                            theta_sep = theta(:,sep_span);
                            f_sep = f(:,sep_span);
                            [~,sep_order] = sort(max(abs(f_sep),[],1),"ascend");

                            x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                            y = [EQUILIBRIUM(2),r_sep(2,sep_order)];
                            z = [0,theta_sep(plotted_outputs(iPlot),sep_order)];
                            plot3(ax{1,iPlot},x,y,z,'.-')
                        end
                        plot3(ax{1,iPlot},EQUILIBRIUM(1),EQUILIBRIUM(2),0,'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end

            end
        end
        %-----------------------------------------------------------------%
        function ax = plot_energy(obj,r,E,sep_id,f)
            EQUILIBRIUM = [0;0];

            modes = obj.Model.reduced_modes;
            num_modes = length(modes);

            if ~exist("r","var")
                r = obj.reduced_displacement;
                E = obj.potential_energy;
                sep_id = obj.static_equilibrium_path_id;
                f = obj.restoring_force;
            end

            num_seps = max(sep_id);
            
            
            switch num_modes
                case 1
                    figure
                    ax = axes(gcf);
                    
                    xlabel("r_" + modes(1))
                    ylabel("V")
                    box on

                    hold on
                    for iSep = 1:num_seps
                        sep_span = sep_id==iSep;
                        r_sep = r(:,sep_span);
                        E_sep = E(:,sep_span);
                        [~,sep_order] = sort(abs(f(1,sep_span)),"ascend");

                        x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                        y = [0,E_sep(1,sep_order)];
                        plot(x,y,'.-')
                    end
                    plot(EQUILIBRIUM,0,'k.',"MarkerSize",8)
                    hold off

                case 2
                    figure
                    ax = axes(gcf);

                    xlabel("r_" + modes(1))
                    ylabel("r_" + modes(2))
                    zlabel("V")
                    box on

                    hold on
                    for iSep = 1:num_seps
                        sep_span = sep_id==iSep;
                        r_sep = r(:,sep_span);
                        E_sep = E(:,sep_span);
                        f_sep = f(:,sep_span);
                        [~,sep_order] = sort(max(abs(f_sep),[],1),"ascend");

                        x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                        y = [EQUILIBRIUM(2),r_sep(2,sep_order)];
                        z = [0,E_sep(1,sep_order)];
                        plot3(x,y,z,'.-')
                    end
                    plot3(EQUILIBRIUM(1),EQUILIBRIUM(2),0,'k.',"MarkerSize",8)
                    hold off
            end
        end
        %-----------------------------------------------------------------%
        function ax = plot_stiffness(obj,plotted_outputs,ax,K_array,r,sep_id)
            EQUILIBRIUM = [0;0];

            num_plots = length(plotted_outputs);
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);


            if nargin == 2
                r = obj.reduced_displacement;
                K_array = obj.tangent_stiffness;
                sep_id = obj.static_equilibrium_path_id;
            end
            K_lin = obj.Model.stiffness;
            sparsity_pattern = K_array.nonzero_indicies(plotted_outputs,:);

            num_seps = max(sep_id);

            if ~exist("ax","var")
                ax = cell(1,num_plots);
                figure
                tiledlayout("flow")
                for iPlot = 1:num_plots

                    ax{1,iPlot} = nexttile;
                end
            end

            if ~iscell(ax)
                ax_temp = ax;
                ax = cell(1,length(ax_temp));
                for iAx = 1:length(ax_temp)
                    ax{1,iAx} = ax_temp(iAx);
                end
                
            end
            

            switch num_modes
                case 1
                    
                    for iPlot = 1:num_plots
                        K_i = K_array.get_component(plotted_outputs(iPlot));
                        row = sparsity_pattern(iPlot,1);
                        col = sparsity_pattern(iPlot,2);

                        xlabel(ax{1,iPlot},"r_" + modes)
                        ylabel(ax{1,iPlot},"K_{(" + row + ',' + col + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [EQUILIBRIUM(1),r(1,sep_span)];
                            y = [K_lin(row,col),K_i(1,sep_span)];
                            plot(ax{1,iPlot},x,y,'.-')
                        end
                        plot(ax{1,iPlot},EQUILIBRIUM(1),K_lin(row,col),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end
                case 2
              

                    for iPlot = 1:num_plots
                        K_i = K_array.get_component(plotted_outputs(iPlot));
                        row = sparsity_pattern(iPlot,1);
                        col = sparsity_pattern(iPlot,2);

                        xlabel(ax{1,iPlot},"r_" + modes(1))
                        ylabel(ax{1,iPlot},"r_" + modes(2))
                        zlabel(ax{1,iPlot},"K_{(" + row + ',' + col + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [EQUILIBRIUM(1),r(1,sep_span)];
                            y = [EQUILIBRIUM(2),r(2,sep_span)];
                            z = [K_lin(row,col),K_i(1,sep_span)];
                            plot3(ax{1,iPlot},x,y,z,'.-')
                        end
                        plot3(ax{1,iPlot},EQUILIBRIUM(1),EQUILIBRIUM(2),K_lin(row,col),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end

            end

        end
        %-----------------------------------------------------------------%
        function ax = plot_perturbation(obj,plotted_outputs,ax)
            EQUILIBRIUM = [0;0];
            

            num_plots = size(plotted_outputs,1);
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);


            f = obj.restoring_force;
            r = obj.reduced_displacement;
            perturbation_disp = obj.perturbation_displacement;
            sep_id = obj.static_equilibrium_path_id;

            %
            stiffness = obj.Model.stiffness;
            mass = obj.Model.mass;
            L_evec = obj.Model.low_frequency_eigenvectors;
            lambda = obj.perturbation_scale_factor;
            perturbation_0 = lambda*(stiffness\(mass*L_evec));
            %
            num_seps = max(sep_id);

            if ~exist("ax","var")
                ax = cell(1,num_plots);
                figure
                tiledlayout("flow")
                for iPlot = 1:num_plots
                    ax{1,iPlot} = nexttile;
                end
            end

            if ~iscell(ax)
                ax_temp = ax;
                ax = cell(1,length(ax_temp));
                for iAx = 1:length(ax_temp)
                    ax{1,iAx} = ax_temp(iAx);
                end
                
            end
            

            switch num_modes
                case 1
                    
                    for iPlot = 1:num_plots
                        plot_index = plotted_outputs(iPlot,:);
                        perturbation_i = squeeze(perturbation_disp(plot_index(1),plot_index(2),:))';
                        perturbation_i_0 = perturbation_0(plot_index(1),plot_index(2));

                        xlabel(ax{1,iPlot},"r_" + modes)
                        ylabel(ax{1,iPlot},"x_{(" + plot_index(1) + ',' + plot_index(2) + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            r_sep = r(:,sep_span);
                            f_sep = f(:,sep_span);
                            perturbation_sep = perturbation_i(:,sep_span);
                            [~,sep_order] = sort(max(abs(f_sep),[],1),"ascend");
                            x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                            y = [perturbation_i_0,perturbation_sep(1,sep_order)];
                            plot(ax{1,iPlot},x,y,'.-')
                        end
                        plot(ax{1,iPlot},EQUILIBRIUM(1),perturbation_i_0,'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end
                case 2
              

                    for iPlot = 1:num_plots
                        plot_index = plotted_outputs(iPlot,:);
                        perturbation_i = squeeze(perturbation_disp(plot_index(1),plot_index(2),:))';


                        xlabel(ax{1,iPlot},"r_" + modes(1))
                        ylabel(ax{1,iPlot},"r_" + modes(2))
                        zlabel(ax{1,iPlot},"G_{(" + plot_index(1) + ',' + plot_index(2) + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [EQUILIBRIUM(1),r(1,sep_span)];
                            y = [EQUILIBRIUM(2),r(2,sep_span)];
                            z = [h_coupling_gradient_0(plot_index(1),plot_index(2)),perturbation_i(1,sep_span)];
                            plot3(ax{1,iPlot},x,y,z,'.-')
                        end
                        plot3(ax{1,iPlot},EQUILIBRIUM(1),EQUILIBRIUM(2),h_coupling_gradient_0(plot_index(1),plot_index(2)),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end

            end
        end
        %-----------------------------------------------------------------%
        function ax = plot_h_stiffness(obj,plotted_outputs,ax,h_stiffness,r,sep_id)
            EQUILIBRIUM = [0;0];

            num_plots = size(plotted_outputs,1);
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);


            if nargin == 2
                r = obj.reduced_displacement;
                h_stiffness = obj.low_frequency_stiffness;
                sep_id = obj.static_equilibrium_path_id;
            end
            h_stiffness_0 = obj.Dynamic_Validation_Data.h_stiffness_0;
            


            num_seps = max(sep_id);

            if ~exist("ax","var")
                ax = cell(1,num_plots);
                figure
                tiledlayout("flow")
                for iPlot = 1:num_plots
                    ax{1,iPlot} = nexttile;
                end
            end

            if ~iscell(ax)
                ax_temp = ax;
                ax = cell(1,length(ax_temp));
                for iAx = 1:length(ax_temp)
                    ax{1,iAx} = ax_temp(iAx);
                end
                
            end
            

            switch num_modes
                case 1
                    
                    for iPlot = 1:num_plots
                        plot_index = plotted_outputs(iPlot,:);
                        h_stiffness_i = squeeze(h_stiffness(plot_index(1),plot_index(2),:))';

                        xlabel(ax{1,iPlot},"r_" + modes)
                        ylabel(ax{1,iPlot},"D_{(" + plot_index(1) + ',' + plot_index(2) + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [EQUILIBRIUM(1),r(1,sep_span)];
                            y = [h_stiffness_0(plot_index(1),plot_index(2)),h_stiffness_i(1,sep_span)];
                            plot(ax{1,iPlot},x,y,'.-')
                        end
                        plot(ax{1,iPlot},EQUILIBRIUM(1),h_stiffness_0(plot_index(1),plot_index(2)),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end
                case 2
              

                    for iPlot = 1:num_plots
                        plot_index = plotted_outputs(iPlot,:);
                        h_stiffness_i = squeeze(h_stiffness(plot_index(1),plot_index(2),:))';


                        xlabel(ax{1,iPlot},"r_" + modes(1))
                        ylabel(ax{1,iPlot},"r_" + modes(2))
                        zlabel(ax{1,iPlot},"D_{(" + plot_index(1) + ',' + plot_index(2) + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [EQUILIBRIUM(1),r(1,sep_span)];
                            y = [EQUILIBRIUM(2),r(2,sep_span)];
                            z = [h_stiffness_0(plot_index(1),plot_index(2)),h_stiffness_i(1,sep_span)];
                            plot3(ax{1,iPlot},x,y,z,'.-')
                        end
                        plot3(ax{1,iPlot},EQUILIBRIUM(1),EQUILIBRIUM(2),h_stiffness_0(plot_index(1),plot_index(2)),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end

            end

        end
        %-----------------------------------------------------------------%
        function ax = plot_h_coupling_gradient(obj,plotted_outputs,ax,h_coupling_gradient,r,sep_id)
            EQUILIBRIUM = [0;0];

            num_plots = size(plotted_outputs,1);
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);


            if nargin == 2
                r = obj.reduced_displacement;
                h_coupling_gradient = obj.low_frequency_coupling_gradient;
                sep_id = obj.static_equilibrium_path_id;
                f = obj.restoring_force;
            end
            h_coupling_gradient_0 = obj.Dynamic_Validation_Data.h_coupling_gradient_0;

            num_seps = max(sep_id);

            if ~exist("ax","var")
                ax = cell(1,num_plots);
                figure
                tiledlayout("flow")
                for iPlot = 1:num_plots
                    ax{1,iPlot} = nexttile;
                end
            end

            if ~iscell(ax)
                ax_temp = ax;
                ax = cell(1,length(ax_temp));
                for iAx = 1:length(ax_temp)
                    ax{1,iAx} = ax_temp(iAx);
                end
                
            end
            

            switch num_modes
                case 1
                    
                    for iPlot = 1:num_plots
                        plot_index = plotted_outputs(iPlot,:);
                        h_coupling_gradient_i = squeeze(h_coupling_gradient(plot_index(1),plot_index(2),:))';

                        xlabel(ax{1,iPlot},"r_" + modes)
                        ylabel(ax{1,iPlot},"G_{(" + plot_index(1) + ',' + plot_index(2) + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            r_sep = r(:,sep_span);
                            f_sep = f(:,sep_span);
                            h_coupling_gradient_sep = h_coupling_gradient_i(:,sep_span);
                            [~,sep_order] = sort(max(abs(f_sep),[],1),"ascend");

                            x = [EQUILIBRIUM(1),r_sep(1,sep_order)];
                            y = [h_coupling_gradient_0(plot_index(1),plot_index(2)),h_coupling_gradient_sep(1,sep_order)];
                            plot(ax{1,iPlot},x,y,'.-')
                        end
                        plot(ax{1,iPlot},EQUILIBRIUM(1),h_coupling_gradient_0(plot_index(1),plot_index(2)),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end
                case 2
              

                    for iPlot = 1:num_plots
                        plot_index = plotted_outputs(iPlot,:);
                        h_coupling_gradient_i = squeeze(h_coupling_gradient(plot_index(1),plot_index(2),:))';


                        xlabel(ax{1,iPlot},"r_" + modes(1))
                        ylabel(ax{1,iPlot},"r_" + modes(2))
                        zlabel(ax{1,iPlot},"G_{(" + plot_index(1) + ',' + plot_index(2) + ")}")
                        box(ax{1,iPlot},"on")

                        hold(ax{1,iPlot},"on")
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [EQUILIBRIUM(1),r(1,sep_span)];
                            y = [EQUILIBRIUM(2),r(2,sep_span)];
                            z = [h_coupling_gradient_0(plot_index(1),plot_index(2)),h_coupling_gradient_i(1,sep_span)];
                            plot3(ax{1,iPlot},x,y,z,'.-')
                        end
                        plot3(ax{1,iPlot},EQUILIBRIUM(1),EQUILIBRIUM(2),h_coupling_gradient_0(plot_index(1),plot_index(2)),'k.',"MarkerSize",8)
                        hold(ax{1,iPlot},"off")
                    end

            end

        end
        %-----------------------------------------------------------------%
        function ax = plot_stress_manifold(obj,f,r,sep_id)
            EQUILIBRIUM = [0;0];
            modes = obj.Model.reduced_modes;
            num_modes = length(modes);

            if ~exist("f","var")
                r = obj.reduced_displacement;
                f = obj.restoring_force;
                sep_id = obj.static_equilibrium_path_id;
            end

            num_seps = max(sep_id);
            
            switch num_modes
                case 1
                    figure
                    ax = axes(gcf);
                    xlabel("f_" + modes(1))
                    ylabel("r_" + modes(1))
                    box on

                    hold on
                    for iSep = 1:num_seps
                        sep_span = sep_id==iSep;
                        
                        x = [0,f(1,sep_span)];
                        y = [EQUILIBRIUM(1),r(1,sep_span)];
                        plot(x,y,'.-')
                    end
                    plot(0,EQUILIBRIUM(1),'k.',"MarkerSize",8)
                    hold off

                case 2
                    figure
                    tiledlayout(1,num_modes)
                    ax = cell(1,num_modes);
                    for iMode = 1:num_modes
                        ax{1,iMode} = nexttile;
                        xlabel("f_" + modes(1))
                        ylabel("f_" + modes(2))
                        zlabel("r_" + modes(iMode))
                        box on

                        hold on
                        for iSep = 1:num_seps
                            sep_span = sep_id==iSep;
                            x = [0,f(1,sep_span)];
                            y = [0,f(2,sep_span)];
                            z = [EQUILIBRIUM(iMode),r(iMode,sep_span)];
                            plot3(x,y,z,'.-')
                        end
                        plot3(0,0,EQUILIBRIUM(iMode),'k.',"MarkerSize",8)
                        hold off
                    end
            end
        end
        %-----------------------------------------------------------------%
    end
end