classdef Dynamic_System
    %Defines the system to be analysed. Can have known or unknowns EoMs
    properties
        system_name
        system_type
        Static_Options
        num_dof
        
        energy_limit
        fitting_energy_limit
        

        Calibration_Options
        calibrated_forces
        calibrated_degree_limits
        
        node_mapping
        mass
        stiffness
        
        reduced_modes
        reduced_eigenvalues
        reduced_eigenvectors 
        
        low_frequency_modes
        low_frequency_eigenvalues
        low_frequency_eigenvectors

        Parameters
    end
    methods
        function obj = Dynamic_System(name,e_lim,modes,Calibration_Opts,Static_Opts,varargin)
            model_init_start = tic;

            %Optional argumanents
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);
            
            load_cache = true; %load prestored mass and stiffness matrices if they exist

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "load_cache"
                        load_cache = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            

            obj.system_name = name;
            obj.energy_limit = e_lim;
            obj.reduced_modes = modes;
            
            
            %----------------------- Environment Setup -------------------%
            if ~isfolder("data\logs")
                mkdir("data\logs")
            end

            log_id = fopen("data\logs\log.txt","w");
            fclose(log_id);
            logger(name + ": " + string(datetime),10)
            logger("",1)

            if isfolder('temp')
                rmdir('temp','s')
            end
            mkdir('temp')
            %-------------------------------------------------------------%
            GEOMETRY_FILE_PATH = "geometry\" + obj.system_name + "\" + obj.system_name;
            if isfile(GEOMETRY_FILE_PATH + ".inp")
                obj.system_type = "indirect";
            elseif isfile(GEOMETRY_FILE_PATH + ".m")
                obj.system_type = "direct";
                Analytic_Eom = load_analytic_system(GEOMETRY_FILE_PATH);
                obj.Parameters = Analytic_Eom.Parameters;
            end

            if obj.system_type == "direct"
                Static_Opts.static_solver = "matlab";
            end
            
            %Update and set optional static solver settings
            obj = obj.update_static_opts(Static_Opts);

            %Update and set optional calibration settings
            obj = obj.update_calibration_opts(Calibration_Opts);
            
            %Find mass and stiffness matricies and find eigenvectors
            matrix_time_start = tic;
            obj = obj.eigenanalysis(load_cache);
            matrix_time = toc(matrix_time_start);
            log_message = sprintf("Eigenvectors: %.1f seconds" ,matrix_time);
            logger(log_message,2)

            %Find single modal forcing required to reach energy limit
            calibration_time_start = tic;
            obj = obj.update_static_opts(Calibration_Opts.Static_Opts);
            obj = obj.calibrate_mode(modes);
            obj = obj.update_static_opts(Static_Opts);
            calibration_time = toc(calibration_time_start);
            log_message = sprintf("Mode Calibration: %.1f seconds" ,calibration_time);
            logger(log_message,2)
            %--------------------%

            model_init_time = toc(model_init_start);
            log_message = sprintf("Model Initialised: %.1f seconds" ,model_init_time);
            logger(log_message,1)
        end
        %-----------------------------------------------------------------%
        function obj = update_static_opts(obj,Static_Opts)
            Default_Static_Opts = read_default_options("static");            
            obj.Static_Options = update_options(Default_Static_Opts,obj.Static_Options,Static_Opts);
        end
        %-----------------------------------------------------------------%
        function obj = update_calibration_opts(obj,Calibration_Opts)
            Default_Calibration_Opts = read_default_options("calibration"); 
            

            New_Calibration_Opts = update_options(Default_Calibration_Opts,obj.Calibration_Options,Calibration_Opts);
            obj.Calibration_Options = New_Calibration_Opts;
            obj.fitting_energy_limit = obj.energy_limit*New_Calibration_Opts.energy_overfit;
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %%% Setup
        %-----------------------------------------------------------------%
        function obj = eigenanalysis(obj,load_cache)
            %extract mass and stiffness matricies and solve generalised
            %eigenvalue problem
            GEOMETRY_PATH = "geometry\" + obj.system_name + "\";
            switch obj.Static_Options.static_solver
                case "abaqus"
                    MATRIX_PATH = GEOMETRY_PATH + "matrices";
                    matrices_loaded = isfile(MATRIX_PATH + ".mat") && load_cache;
                    if matrices_loaded
                        load(MATRIX_PATH,"M","K","node_map")

                        logger("Matrices Loaded",3)
                    else
                        [M,K,node_map] = matrices_abaqus(obj.system_name);
                        save(MATRIX_PATH,"M","K","node_map")
                    end

                case "matlab"
                    Analytic_Eom = load_analytic_system(GEOMETRY_PATH + obj.system_name);
                    M = Analytic_Eom.linear_mass;
                    K = Analytic_Eom.linear_stiffness;

                    dofs = length(M);
                    node_map = [(1:dofs)',(1:dofs)'];

            end
            obj.num_dof = size(node_map,1);
            obj.node_mapping = node_map;
            obj.mass = M;
            obj.stiffness = K;
            
            
            switch obj.Static_Options.additional_data
                case "perturbation"
                    r_modes = obj.reduced_modes;
                    L_modes = 1:obj.Static_Options.num_validation_modes;
                    L_modes(ismember(L_modes,r_modes)) = [];
                    obj.low_frequency_modes = L_modes;
                    h_modes = [r_modes,L_modes];

                    [eVec,eVal] = eigs(K,M,max(h_modes),"smallestabs");

                    eVal_L = eVal(L_modes,L_modes)*ones(length(L_modes),1);
                    eVec_L = eVec(:,L_modes);

                    obj.low_frequency_eigenvalues = eVal_L;
                    obj.low_frequency_eigenvectors = eVec_L;

                otherwise
                    r_modes = obj.reduced_modes;
                    [eVec,eVal] = eigs(K,M,max(r_modes),"smallestabs");

            end

            eVal_r = eVal(r_modes,r_modes)*ones(length(r_modes),1);
            eVec_r = eVec(:,r_modes);

            obj.reduced_eigenvalues = eVal_r;
            obj.reduced_eigenvectors = eVec_r;
        end
        %-----------------------------------------------------------------%
        function obj = calibrate_mode(obj,modes)
            %finds static forces that reach the potential energy limit
            Calibration_Opts = obj.Calibration_Options;

            %check if already calibrated
            GEOMETRY_PATH = "geometry\" + obj.system_name + "\";
            if isfile(GEOMETRY_PATH + "force_calibration.mat")
                load(GEOMETRY_PATH + "force_calibration.mat","Force_Calibration");
            end

            if isfile(GEOMETRY_PATH + "force_calibration.mat") && isequal(Force_Calibration.Parameters,obj.Parameters)
                calibrated_energy = Force_Calibration.energy_limit;
                if ismember(obj.energy_limit,calibrated_energy)
                    calibration_id = find(calibrated_energy == obj.energy_limit);
                    calibrated_modes = Force_Calibration.calibrated_modes{1,calibration_id};
                    uncalibrated_modes = setdiff(modes,calibrated_modes);
                else
                    calibration_id = length(calibrated_energy) + 1;
                    Force_Calibration.energy_limit(calibration_id) = obj.energy_limit;
                    uncalibrated_modes = modes;
                    calibrated_modes = [];
                end
            else
                uncalibrated_modes = modes;
                calibrated_modes = [];
                Force_Calibration.energy_limit = obj.energy_limit;
                Force_Calibration.force_limit = {};
                Force_Calibration.calibrated_modes = {};
                calibration_id = 1;
            end
            
            %calibrate force scale factors
            num_calibrated_modes = length(calibrated_modes);
            num_uncalibrated_modes = length(uncalibrated_modes);
           
            
            r_modes = obj.reduced_modes;
            num_r_modes = length(r_modes);
            r_eigenvalues = obj.reduced_eigenvalues;

            num_matching_calibrated_modes = length(intersect(r_modes,calibrated_modes));
            
            log_message = sprintf("%u/%u modes precalibrated",[num_matching_calibrated_modes,num_r_modes]);
            logger(log_message,3)

            initial_force_ratio = zeros(num_r_modes,num_uncalibrated_modes*2);
            for iMode = 1:num_uncalibrated_modes
                mode = uncalibrated_modes(iMode);
                mode_index = mode == modes;

                force_ratio = zeros(num_r_modes,2);
                force_ratio(mode_index,:) = [1,-1];

                %start with linear approximation
                force_scale_factor = Calibration_Opts.calibration_scale_factor*sqrt(2*r_eigenvalues(mode_index)*obj.fitting_energy_limit);
                initial_force_ratio(:,[2*iMode-1,2*iMode]) = force_ratio*force_scale_factor;

            end

            if num_uncalibrated_modes > 0
                [r,x,f,E,sep_id] = obj.add_sep(initial_force_ratio);
            end

            for iMode = 1:num_uncalibrated_modes
                mode = uncalibrated_modes(iMode);
                mode_index = mode == modes;

                num_seps = 2;
                r_limit = zeros(1,num_seps);
                f_limit = zeros(1,num_seps);
                for iSep = 1:num_seps
                    sep_span = sep_id == (iSep+2*(iMode-1));
                    E_sep = E(sep_span);
                    r_sep = r(mode_index,sep_span);
                    f_sep = f(mode_index,sep_span);

                    E_diff = E_sep - obj.energy_limit;
                    E_upper = E_diff;
                    E_lower = E_diff;

                    E_upper(E_upper < 0) = inf;
                    E_lower(E_lower > 0) = -inf;

                    [~,min_index] = max(E_lower);
                    [~,max_index] = min(E_upper);
                    bound_index = [min_index,max_index];

                    r_limit(1,iSep) = interp1(E_sep(bound_index),r_sep(bound_index),obj.energy_limit);
                    f_limit(1,iSep) = interp1(E_sep(bound_index),f_sep(bound_index),obj.energy_limit);
                end
                current_seps = [1,2]+2*(iMode-1);
                % Min_Degree_Data = find_degree_limits(obj,r,x,f,E,sep_id,current_seps);


                Force_Calibration.force_limit{1,calibration_id}(iMode+num_calibrated_modes,:) = f_limit;
                Force_Calibration.calibrated_modes{1,calibration_id}(iMode+num_calibrated_modes,:) = mode;
                % Force_Calibration.min_degree_data{1,calibration_id}{iMode+num_calibrated_modes} = Min_Degree_Data;
                % f_limit = interp1(E,r,obj.energy_limit);

                % r_limits = [min(r),max(r)];
                %---------------------------------------------------------%
                figure
                tiledlayout(1,2)
                nexttile
                box on
                xlabel("r_" + mode)
                ylabel("f_" + mode)

                hold on
                for iSep = 1:2
                    sep_span = sep_id == (iSep+2*(iMode-1));
                    plot([0,r(mode_index,sep_span)],[0,f(mode_index,sep_span)],'.-')
                end

                for iSep = 1:2
                    plot(gca().XLim,f_limit(iSep)*[1,1],'k-')
                end
                plot(0,0,'k.','MarkerSize',10)
                hold off

                nexttile
                box on
                xlabel("r_" + mode)
                ylabel("V")
                
                hold on
                
                for iSep = 1:2
                    sep_span = sep_id == (iSep+2*(iMode-1));
                    plot([0,r(mode_index,sep_span)],[0,E(sep_span)],'.-')
                end
                plot(gca().XLim,obj.fitting_energy_limit*[1,1],'k-')
                plot(0,0,'k.','MarkerSize',10)
                hold off
                %---------------------------------------------------------%
            end
            Force_Calibration.Parameters = obj.Parameters;
            save(GEOMETRY_PATH + "force_calibration","Force_Calibration")
            
            calibrated_modes = Force_Calibration.calibrated_modes{1,calibration_id};
            obj.calibrated_forces = zeros(num_r_modes,2);
            % obj.calibrated_degree_limits = Force_Calibration.min_degree_data{1,calibration_id};

            for iMode = 1:num_r_modes
                mode = r_modes(iMode);
                obj.calibrated_forces(iMode,:) = Force_Calibration.force_limit{1,calibration_id}(mode == calibrated_modes,:)*Calibration_Opts.force_overcalibration;
                % obj.calibrated_degree_limits{iMode}.force_applied_force = obj.calibrated_degree_limits{iMode}.force_applied_force./obj.calibrated_forces(iMode,:)';
                % obj.calibrated_degree_limits{iMode}.disp_applied_force = obj.calibrated_degree_limits{iMode}.disp_applied_force./obj.calibrated_forces(iMode,:)';
            end
            
            

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %%% Simulation
        %-----------------------------------------------------------------%
        function [reduced_disp,condensed_disp,restoring_force,energy,sep_id,additional_data] = ...
                add_sep(obj,force_ratio,additional_data_type,clean_data)
            % Add the Static Equilibrium Path (SEP) corresponding to the force ratio
            % Force ratio also defines initial magnitude and direction
            if ~exist("additional_data_type","var")
                additional_data_type = "none";
            end

            if ~exist("clean_data","var")
                clean_data = 0;
            end

            Static_Opts = obj.Static_Options;
            num_loadcases = Static_Opts.num_loadcases;


            switch Static_Opts.static_solver
                case "abaqus"
                    max_inc = Static_Opts.maximum_step_increments;
                    reset_temp_directory()
                    if Static_Opts.max_parallel_jobs > 1
                        abaqus_start = tic;
                        max_parallel_jobs = Static_Opts.max_parallel_jobs;
                        mininum_job_loadcases = Static_Opts.minimum_job_loadcases;
                        force_ratio_groups = split_abaqus_jobs(force_ratio,Static_Opts.num_loadcases,max_parallel_jobs,mininum_job_loadcases);
                        num_parallel_jobs = size(force_ratio_groups,2);

                        log_message = sprintf("%i SEPs over %i jobs" ,[size(force_ratio,2),num_parallel_jobs]);
                        logger(log_message,3)
                        logger("---",3)
                        
                        reduced_disp_cell = cell(1,num_parallel_jobs);
                        condensed_disp_cell = cell(1,num_parallel_jobs);
                        restoring_force_cell = cell(1,num_parallel_jobs);
                        energy_cell = cell(1,num_parallel_jobs);
                        additional_data_cell = cell(1,num_parallel_jobs);
                        sep_id_cell = cell(1,num_parallel_jobs);

                        parfor iJob = 1:num_parallel_jobs
                            job_force = force_ratio_groups{1,iJob};
                            [job_r,job_theta,job_f,job_E,job_additional_data,job_sep_id] = ...
                                add_sep_abaqus(job_force,num_loadcases,Static_Opts,max_inc,additional_data_type,clean_data,obj,iJob);
                            
                            reduced_disp_cell{1,iJob} = job_r;
                            condensed_disp_cell{1,iJob} = job_theta;
                            restoring_force_cell{1,iJob} = job_f;
                            energy_cell{1,iJob} = job_E;
                            additional_data_cell{1,iJob} = job_additional_data; 
                            sep_id_cell{1,iJob} = job_sep_id;
                        end
                        
                        reduced_disp = [reduced_disp_cell{1,:}];
                        condensed_disp = [condensed_disp_cell{1,:}];
                        restoring_force = [restoring_force_cell{1,:}];
                        energy = [energy_cell{1,:}];
                        additional_data = cat(3,additional_data_cell{1,:});
                        
                        
                        
                        sep_id_shift = max(sep_id_cell{1,1});
                        for iSep = 2:num_parallel_jobs
                            sep_id_cell{1,iSep} = sep_id_cell{1,iSep} + sep_id_shift;
                            sep_id_shift = max(sep_id_cell{1,iSep});
                        end
                        sep_id = [sep_id_cell{1,:}];
                        
                        abaqus_time = toc(abaqus_start);
                        logger("---",3)
                        log_message = sprintf("Total FE time: %.1f seconds" ,abaqus_time);
                        logger(log_message,3)
                        
                    else
                        [reduced_disp,condensed_disp,restoring_force,energy,additional_data,sep_id] = ...
                        add_sep_abaqus(force_ratio,num_loadcases,Static_Opts,max_inc,additional_data_type,clean_data,obj,1);
                    end


                   
                case "matlab"
                    [reduced_disp,condensed_disp,restoring_force,energy,additional_data,sep_id] = ...
                        add_sep_matlab(force_ratio,num_loadcases,additional_data_type,obj);
            end
        end
        %-----------------------------------------------------------------%
        function [reduced_disp,condensed_disp,restoring_force,energy,additional_data] = ...
                add_point(obj,applied_force,additional_data_type)
            
            if ~exist("additional_data_type","var")
                additional_data_type = "none";
            end


            Static_Opts = obj.Static_Options;
            switch Static_Opts.static_solver
                case "abaqus"
                    
                    reset_temp_directory()
                    max_inc = Static_Opts.maximum_step_increments*Static_Opts.num_loadcases;

                    if Static_Opts.max_parallel_jobs > 1
                        
                        abaqus_start = tic;
                        max_parallel_jobs = Static_Opts.max_parallel_jobs;
                        mininum_job_loadcases = Static_Opts.minimum_job_loadcases;
                        force_groups = split_abaqus_jobs(applied_force,1,max_parallel_jobs,mininum_job_loadcases);
                        num_parallel_jobs = size(force_groups,2);

                        log_message = sprintf("%i loadcases over %i jobs" ,[size(applied_force,2),num_parallel_jobs]);
                        logger(log_message,3)
                        logger("---",3)

                        reduced_disp_cell = cell(1,num_parallel_jobs);
                        condensed_disp_cell = cell(1,num_parallel_jobs);
                        restoring_force_cell = cell(1,num_parallel_jobs);
                        energy_cell = cell(1,num_parallel_jobs);
                        additional_data_cell = cell(1,num_parallel_jobs);

                        parfor iJob = 1:num_parallel_jobs
                            job_force = force_groups{1,iJob};
                            [job_r,job_theta,job_f,job_E,job_additional_data] = ...
                                add_point_abaqus(job_force,max_inc,additional_data_type,obj,iJob);

                            reduced_disp_cell{1,iJob} = job_r;
                            condensed_disp_cell{1,iJob} = job_theta;
                            restoring_force_cell{1,iJob} = job_f;
                            energy_cell{1,iJob} = job_E;
                            additional_data_cell{1,iJob} = job_additional_data;
 
                        end

                        reduced_disp = [reduced_disp_cell{1,:}];
                        condensed_disp = [condensed_disp_cell{1,:}];
                        restoring_force = [restoring_force_cell{1,:}];
                        energy = [energy_cell{1,:}];
                        additional_data = cat(3,additional_data_cell{1,:});
                        
                        abaqus_time = toc(abaqus_start);
                        logger("---",3)
                        log_message = sprintf("Total FE time: %.1f seconds" ,abaqus_time);
                        logger(log_message,3)
                    else

                        [reduced_disp,condensed_disp,restoring_force,energy,additional_data] = ...
                            add_point_abaqus(applied_force,max_inc,additional_data_type,obj,1);
                    end
                case "matlab"
                    [reduced_disp,condensed_disp,restoring_force,energy,additional_data] = ...
                        add_point_matlab(applied_force,additional_data_type,obj);
            end
        end
        %-----------------------------------------------------------------%
        function [t,x,x_dot,energy] = dynamic_simulation(obj,x_0,x_dot_0,f_r_0,period,num_periods,min_incs,initial_time,FE_Force_Data,job_id)
            if ~exist("initial_time","var")
                initial_time = zeros(1,size(x_0,2));
            end
            if ~exist("FE_Force_Data","var")
                FE_Force_Data = [];
            end

            Static_Opts = obj.Static_Options;
            if ~exist("job_id","var")
                job_id = [];
                if Static_Opts.static_solver == "abaqus"
                    reset_temp_directory()
                end
            end
            
            switch Static_Opts.static_solver
                case "abaqus"
                    

                    if Static_Opts.max_parallel_jobs > 1
                        abaqus_start = tic;
                        num_parallel_jobs = size(f_r_0,2);

                        log_message = sprintf("%i dynamic simulations over %i jobs" ,[size(f_r_0,2),num_parallel_jobs]);
                        logger(log_message,3)
                        logger("---",3)

                        t = cell(1,num_parallel_jobs);
                        x = cell(1,num_parallel_jobs);
                        x_dot = cell(1,num_parallel_jobs);
                        energy = cell(1,num_parallel_jobs);
                        parfor (iJob = 1:num_parallel_jobs,Static_Opts.max_parallel_jobs)
                        % for iJob = 1:num_parallel_jobs
                            [t_job,x_job,x_dot_job,energy_job] = dynamic_simulation_abaqus(x_0(:,iJob),x_dot_0(:,iJob),f_r_0(:,iJob),period(iJob),num_periods,min_incs(iJob),initial_time(iJob),FE_Force_Data,obj,iJob);

                            t{1,iJob} = t_job;
                            x{1,iJob} = x_job;
                            x_dot{1,iJob} = x_dot_job;
                            energy{1,iJob} = energy_job;
                        end

                        abaqus_time = toc(abaqus_start);
                        logger("---",3)
                        log_message = sprintf("Total FE time: %.1f seconds" ,abaqus_time);
                        logger(log_message,3)
                    else
                        [t,x,x_dot,energy] = dynamic_simulation_abaqus(x_0,x_dot_0,f_r_0,period,num_periods,min_incs,initial_time,FE_Force_Data,obj,job_id);
                    end
                case "matlab"

            end
        end
        %-----------------------------------------------------------------%
        function [stress,stress_labels,section_points] = stress_simulation(obj,x_0,x_dot_0,f_r_0)
            
            Static_Opts = obj.Static_Options;
            % Static_Opts.max_parallel_jobs = 1;
            switch Static_Opts.static_solver
                case "abaqus"
                    reset_temp_directory()

                    if Static_Opts.max_parallel_jobs > 1
                        abaqus_start = tic;
                        num_parallel_jobs = size(f_r_0,2);

                        log_message = sprintf("%i stress simulations over %i jobs" ,[size(f_r_0,2),num_parallel_jobs]);
                        logger(log_message,3)
                        logger("---",3)

                        stress = cell(1,num_parallel_jobs);
                        section_points = cell(1,num_parallel_jobs);
                        stress_labels = cell(1,num_parallel_jobs);
                        parfor (iJob = 1:num_parallel_jobs,Static_Opts.max_parallel_jobs)
                            [stress_job,stress_labels_job,section_points_job] = stress_simulation_abaqus(x_0(:,iJob),x_dot_0(:,iJob),f_r_0(:,iJob),obj,iJob);
                            stress{1,iJob} = stress_job;
                            stress_labels{1,iJob} = stress_labels_job;
                            section_points{1,iJob} = section_points_job;
                        end
                        
                        abaqus_time = toc(abaqus_start);
                        logger("---",3)
                        log_message = sprintf("Total FE time: %.1f seconds" ,abaqus_time);
                        logger(log_message,3)
                    else
                        [stress,stress_labels,section_points] = stress_simulation_abaqus(x_0,x_dot_0,f_r_0,obj,1);
                    end
                case "matlab"
                    

            end
        end
        %-----------------------------------------------------------------%
    
        %-----------------------------------------------------------------%
        %%% Helpers %%%
        %-----------------------------------------------------------------%
        function geometry = load_geometry(obj)
            switch obj.system_type
                case "indirect"
                    G_ID =fopen("geometry\" + obj.system_name+ "\" + obj.system_name + ".inp");
                    geometry = textscan(G_ID,'%s','delimiter','\n');
                    fclose(G_ID);
                    geometry = geometry{1,1};
                case "direct"

            end
        end
        %-----------------------------------------------------------------%
        function obj = get_modal_subset(obj,modes)
            mode_map = ismember(obj.reduced_modes,modes);
            if all(mode_map)
                return
            end

            obj.reduced_modes = modes;
            obj.reduced_eigenvalues = obj.reduced_eigenvalues(mode_map,1);
            obj.reduced_eigenvectors = obj.reduced_eigenvectors(:,mode_map);
        end
    end
end