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

        dof_boundary_conditions
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
        function obj = Dynamic_System(name,e_lim,modes,varargin)
            model_init_start = tic;

            %Optional argumanents
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Calibration_Opts = struct([]);
            Static_Opts = struct([]);
            load_cache = true; %load prestored mass and stiffness matrices if they exist

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "load_cache"
                        load_cache = keyword_values{arg_counter};
                    case "static_opts"
                        Static_Opts = keyword_values{arg_counter};
                    case "calibration_opts"
                        Calibration_Opts = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end


            obj.system_name = name;
            obj.energy_limit = e_lim;



            %----------------------- Environment Setup -------------------%
            if ~isfolder("data\logs")
                mkdir("data\logs")
            end

            log_id = fopen("data\logs\log.txt","w");
            fclose(log_id);
            logger(name + ": " + string(datetime),10)
            logger("",1)

            if isfolder('temp')
                try
                    rmdir('temp','s')
                catch
                    pause(1)
                    rmdir('temp','s')
                end
            end
            mkdir('temp')
            %-------------------------------------------------------------%
            geometry_file_path = "geometry\" + obj.system_name + "\" + obj.system_name;
            if isfile(geometry_file_path + ".inp")
                obj.system_type = "indirect";
            elseif isfile(geometry_file_path + ".m")
                obj.system_type = "direct";
                Analytic_Eom = load_analytic_system(geometry_file_path);
                obj.Parameters = Analytic_Eom.Parameters;
            end

            if obj.system_type == "direct"
                Static_Opts.static_solver = "matlab";
            end


           


            obj.reduced_modes = modes;
            if isfolder(obj.get_data_path)
                rmdir(obj.get_data_path,"s")
            end

            %Update and set optional static solver settings
            obj = obj.update_static_opts(Static_Opts);
            

            %Update and set optional calibration settings
            obj = obj.update_calibration_opts(Calibration_Opts);

            %Extract  data from input file
            obj = obj.set_problem_data;

            %Find mass and stiffness matricies and find eigenvectors
            matrix_time_start = tic;
            obj = obj.eigenanalysis(load_cache);
            matrix_time = toc(matrix_time_start);
            log_message = sprintf("Eigenvectors: %.1f seconds" ,matrix_time);
            logger(log_message,2)
            if obj.energy_limit > 0
                %Find single modal forcing required to reach energy limit
                calibration_time_start = tic;
                obj = obj.update_static_opts(obj.Calibration_Options.Static_Opts);
                obj = obj.calibrate_mode(modes);
                obj.Static_Options = struct([]);
                obj = obj.update_static_opts(Static_Opts);
                calibration_time = toc(calibration_time_start);
                log_message = sprintf("Mode Calibration: %.1f seconds" ,calibration_time);
                logger(log_message,2)
                %--------------------%
            end

          

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
            if ~isstruct(New_Calibration_Opts.Static_Opts)
                New_Calibration_Opts.Static_Opts = struct([]);
            end
            obj.Calibration_Options = New_Calibration_Opts;
            obj.fitting_energy_limit = obj.energy_limit*New_Calibration_Opts.energy_overfit;
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        %%% Setup
        %-----------------------------------------------------------------%
        function obj = set_problem_data(obj)
            if obj.system_type == "direct"
                return
            end
            if obj.Static_Options.static_solver ~= "abaqus"
                return
            end


            geometry = load_geometry(obj);
            mesh_data = get_mesh_data(geometry);

            %------------ save data
            data_path = "geometry\" + obj.system_name + "\mesh_data";
            save(data_path,"mesh_data")
        end
        %-----------------------------------------------------------------%
        function obj = eigenanalysis(obj,load_cache)
            %extract mass and stiffness matricies and solve generalised
            %eigenvalue problem
            obj = model_eigenanalysis(obj,load_cache);
            

        end
        %-----------------------------------------------------------------%
        function obj = calibrate_mode(obj,modes)
            %finds static forces that reach the potential energy limit
            
            
            if obj.Calibration_Options.disable_calibration
                return
            end

            %static settings
            Static_Opts = obj.Static_Options;
            if isstring(Static_Opts.num_loadcases) && Static_Opts.num_loadcases == "auto"
                Static_Opts.num_loadcases = 5;
                obj = obj.update_static_opts(Static_Opts);
            end

            obj = modal_calibration(obj,modes);
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
                    num_dofs = obj.num_dof;
                    num_seps = size(force_ratio,2);
                    [num_sep_input_lines,additional_job_input_lines] = estimate_input_lines(num_dofs,Static_Opts);
                    max_parallel_jobs = maximum_viable_jobs(num_seps,num_sep_input_lines,additional_job_input_lines,Static_Opts);

                    current_pool = gcp("nocreate");
                    if isempty(current_pool) && max_parallel_jobs > 1
                        current_pool = parpool(max_parallel_jobs);
                    end
                    if ~isempty(current_pool) && current_pool.NumWorkers ~= max_parallel_jobs
                        delete(current_pool)
                        if max_parallel_jobs > 1
                            current_pool = parpool(max_parallel_jobs);
                        end
                    end

                    if max_parallel_jobs > 1
                        abaqus_start = tic;
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
                        switch additional_data_type
                            case "stiffness"
                                additional_data = cat(3,additional_data_cell{1,:});
                                % additional_data = additional_data_cell{1,1};
                                % for iJob = 2:num_parallel_jobs
                                %     additional_data = cat(3,additional_data,additional_data_cell{1,iJob});
                                % end
                            otherwise
                                additional_data = cat(3,additional_data_cell{1,:});
                        end



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

                    elseif max_parallel_jobs == 1
                        [reduced_disp,condensed_disp,restoring_force,energy,additional_data,sep_id] = ...
                            add_sep_abaqus(force_ratio,num_loadcases,Static_Opts,max_inc,additional_data_type,clean_data,obj,1);
                    elseif  max_parallel_jobs < 1
                        % go sep by SEP. For larger systems may have to
                        % restart SEPs mid way
                        abaqus_start = tic;

                        reduced_disp_cell = cell(1,num_seps);
                        condensed_disp_cell = cell(1,num_seps);
                        restoring_force_cell = cell(1,num_seps);
                        energy_cell = cell(1,num_seps);
                        additional_data_cell = cell(1,num_seps);
                        sep_id_cell = cell(1,num_seps);

                        for iSep = 1:num_seps
                            [job_r,job_x,job_f,job_E,job_additional_data,job_sep_id] = ...
                                add_sep_abaqus(force_ratio(:,iSep),num_loadcases,Static_Opts,max_inc,additional_data_type,clean_data,obj,iSep);

                            reduced_disp_cell{1,iSep} = job_r;
                            condensed_disp_cell{1,iSep} = job_x;
                            restoring_force_cell{1,iSep} = job_f;
                            energy_cell{1,iSep} = job_E;
                            additional_data_cell{1,iSep} = job_additional_data;
                            sep_id_cell{1,iSep} = job_sep_id;
                        end

                        reduced_disp = [reduced_disp_cell{1,:}];
                        condensed_disp = [condensed_disp_cell{1,:}];
                        restoring_force = [restoring_force_cell{1,:}];
                        energy = [energy_cell{1,:}];
                        switch additional_data_type
                            case "stiffness"
                                additional_data = additional_data_cell{1,1};
                                for iJob = 2:num_seps
                                    additional_data = cat(3,additional_data,additional_data_cell{1,iJob});
                                end
                            otherwise
                                additional_data = cat(3,additional_data_cell{1,:});
                        end



                        sep_id_shift = max(sep_id_cell{1,1});
                        for iSep = 2:num_seps
                            sep_id_cell{1,iSep} = sep_id_cell{1,iSep} + sep_id_shift;
                            sep_id_shift = max(sep_id_cell{1,iSep});
                        end
                        sep_id = [sep_id_cell{1,:}];


                        abaqus_time = toc(abaqus_start);
                        logger("---",3)
                        log_message = sprintf("Total FE time: %.1f seconds" ,abaqus_time);
                        logger(log_message,3)
                    end



                case "matlab"
                    switch Static_Opts.solver_algorithm
                        case "standard"
                            [reduced_disp,condensed_disp,restoring_force,energy,additional_data,sep_id] = ...
                                add_sep_matlab(force_ratio,num_loadcases,additional_data_type,obj);
                        case "riks"
                            [reduced_disp,condensed_disp,restoring_force,energy,additional_data,sep_id] = ...
                                add_sep_matlab_riks(force_ratio,num_loadcases,additional_data_type,obj);
                    end
            end
        end
        %-----------------------------------------------------------------%
        function [reduced_disp,condensed_disp,restoring_force,energy,additional_data] = ...
                add_point(obj,applied_force,additional_data_type,Closest_Point)

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
                        [force_groups,job_index] = split_abaqus_jobs(applied_force,1,max_parallel_jobs,mininum_job_loadcases);
                        num_parallel_jobs = size(force_groups,2);
                        closest_point_group = cell(size(force_groups));
                        if ~isempty(Closest_Point)
                            num_jobs = size(force_groups,2);
                            for iJob = 1:num_jobs
                                Job_Closest_Point.initial_disp = Closest_Point.initial_disp(:,job_index{iJob});
                                Job_Closest_Point.initial_force = Closest_Point.initial_force(:,job_index{iJob});
                                closest_point_group{iJob} = Job_Closest_Point;
                            end
                        end

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
                            Job_Closest_Point = closest_point_group{1,iJob};
                            [job_r,job_theta,job_f,job_E,job_additional_data] = ...
                                add_point_abaqus(job_force,max_inc,additional_data_type,obj,iJob,Job_Closest_Point);

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
                        switch additional_data_type
                            case "stiffness"
                                additional_data = additional_data_cell{1,1};
                                for iJob = 2:num_parallel_jobs
                                    additional_data = cat(3,additional_data,additional_data_cell{1,iJob});
                                end
                            otherwise
                                additional_data = cat(3,additional_data_cell{1,:});
                        end

                        abaqus_time = toc(abaqus_start);
                        logger("---",3)
                        log_message = sprintf("Total FE time: %.1f seconds" ,abaqus_time);
                        logger(log_message,3)
                    else

                        [reduced_disp,condensed_disp,restoring_force,energy,additional_data] = ...
                            add_point_abaqus(applied_force,max_inc,additional_data_type,obj,1,Closest_Point);
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
        function [eom,eom_dz] = get_equation_of_motion(obj,varargin)
            Analytic_Eom = load_analytic_system("geometry\" + obj.system_name+ "\" + obj.system_name);
            %-------------------------------------------------------------------------%
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Damping_Data = [];
            Force_Data = [];

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "damping"
                        Damping_Data = keyword_values{arg_counter};
                    case "forcing"
                        Force_Data = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------%
            conservative = (isempty(Damping_Data) && isempty(Force_Data));

            if conservative
                Eom_Input = Analytic_Eom.get_solver_inputs("free");

                eom = @(t,z) direct_eom(0,z,0,Eom_Input.modal_restoring_force);
                eom_dz = @(t,z) direct_eom_dx(0,z,0,Eom_Input.modal_stiffness);
            else
                switch Damping_Data.damping_type
                    case "rayleigh"
                        damping = get_rayleigh_damping_matrix(Damping_Data,obj);
                        % evec = obj.reduced_eigenvectors;
                        % damping = evec'*damping*evec;
                        %need to transform to modal coordinates
                end
                if ~isempty(Force_Data)
                    error("forced full-order simulation not implemented")
                end
                Eom_Input = Analytic_Eom.get_solver_inputs("free");
                %implement properly later

                eom = @(t,z) direct_forced_eom(t,z,Eom_Input.modal_restoring_force,damping);

            end

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        %%% Helpers %%%
        %-----------------------------------------------------------------%
        function data_path = get_data_path(obj)
            r_modes = obj.reduced_modes;
            mode_id = join(string(r_modes),"");
            data_path = "data\" + obj.system_name + "_" + mode_id + "\model_data\";
        end
        %-------------------------------------------
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