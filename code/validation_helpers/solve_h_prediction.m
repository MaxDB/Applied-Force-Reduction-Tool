function Validation_Sol = solve_h_prediction(Validation_Sol,Solution,Validation_Rom, Validated_BB_Settings)
Validation_Opts = Validation_Sol.Validation_Options;
solution_num = Validated_BB_Settings.solution_num;


num_r_modes = length(Validation_Rom.Model.reduced_modes);
disp_span = 1:num_r_modes;
vel_span = disp_span + num_r_modes;


Solution_Type = Solution.Solution_Type;
orbit_type = Solution_Type.orbit_type;
solution_name = Validation_Rom.data_path + "dynamic_sol_" + solution_num;



%%% Set up h-problem

switch orbit_type
    case "free"
        Eom_Input = Validation_Rom.get_solver_inputs("coco_backbone");
        reduced_eom = @(t,z,zeta) coco_eom(t,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
        %
        Validation_Input = Validation_Rom.get_solver_inputs("h_prediction");
        h_terms = @(r,r_dot,r_ddot) get_h_error_terms(r,r_dot,r_ddot,Validation_Input);

        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("h_analysis",Validated_BB_Settings.Additional_Output);
    case "forced"
        
        Nonconservative_Input = Solution.get_nonconservative_input(Validation_Rom.Model);
        amp = Nonconservative_Input.amplitude;
        Eom_Input = Validation_Rom.get_solver_inputs("coco_frf",Nonconservative_Input);
        reduced_eom = @(t,z,T) coco_forced_eom(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);
        
        Validation_Input = Validation_Rom.get_solver_inputs("forced_h_prediction",Nonconservative_Input);
        h_terms = @(t,r,r_dot,r_ddot,period) get_forced_h_error_terms(t,r,r_dot,r_ddot,amp,period,Validation_Input);

        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("forced_h_analysis",Nonconservative_Input);
end


switch Validation_Opts.validation_algorithm
    case "h_time"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_time_solution(h_terms,t0,omega,num_harmonics);
    case "h_frequency"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_harmonic_balance(h_terms,t0,omega,num_harmonics);
    case "h_infinite_determinant"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_infinite_determinant(h_terms,t0,omega,num_harmonics);
end


orbit_labels = Solution.orbit_labels;
frequency = Solution.frequency;
num_periodic_orbits = length(orbit_labels);

num_jobs = 4;
orbit_groups = split_orbit_jobs(num_periodic_orbits,num_jobs);


Validation_Sols = repelem(Validation_Sol,1,num_jobs);


initial_harmonic = Validation_Opts.initial_harmonic;
Force_Polynomial = Validation_Rom.Force_Polynomial;
time_ranges = zeros(2,num_jobs);
job_time = zeros(1,num_jobs);
num_job_orbits = zeros(1,num_jobs);

parfor iJob = 1:num_jobs
    time_range = [inf,0];
    num_harmonics = initial_harmonic;
    orbit_group = orbit_groups(iJob,:);
    periodic_orbits_jobs = orbit_group(1):orbit_group(2);
    num_periodic_orbits_jobs = size(periodic_orbits_jobs,2);
    num_job_orbits(iJob) = num_periodic_orbits_jobs;    

    job_orbit_labels = orbit_labels(periodic_orbits_jobs);
    job_frequency = frequency(periodic_orbits_jobs);
    for iOrbit = 1:num_periodic_orbits_jobs
        read_data_start = tic;
        job_orbit = periodic_orbits_jobs(iOrbit);
        sol = po_read_solution('',convertStringsToChars(solution_name),job_orbit_labels(iOrbit));
        read_data_time = toc(read_data_start);

        set_up_h_start = tic;
        t0 = sol.tbp';
        x = sol.xbp';
        r = x(disp_span,:);
        r_dot = x(vel_span,:);
        omega = job_frequency(iOrbit);

        switch orbit_type
            case "free"
                x_dot = reduced_eom(t0,x,zeros(size(t0)));
                r_ddot  = x_dot(vel_span,:);
                [h_inertia,h_conv,h_stiff,h_force] = h_terms(r,r_dot,r_ddot); %lots of scope to speed up
            case "forced"

                period = 2*pi/omega;
                x_dot = reduced_eom(t0,x,period);
                r_ddot  = x_dot(vel_span,:);
                [h_inertia,h_conv,h_stiff,h_force] = h_terms(t0,r,r_dot,r_ddot,period);
            otherwise
                h_inertia = [];
                h_conv = [];
                h_stiff = [];
                h_force = [];
                error("")
        end
        %h_inertia * h_ddot  +  h_conv * h_dot  +  h_stiff * h  =  h_force
        set_up_h_time = toc(set_up_h_start);

        validation_eq_terms = {h_inertia,h_conv,h_stiff,h_force};
        r_force = Force_Polynomial.evaluate_polynomial(r);
        solution_converged = 0;
        solve_h_time = 0;
        Validation_Orbit = [];
        while ~solution_converged
            solve_h_start = tic;
            h_frequency = h_solver(validation_eq_terms,t0,omega,num_harmonics);
            [solution_converged,num_harmonics,Validation_Orbit] = check_h_convergence(validation_eq_terms,r_force,h_frequency,t0,omega,num_harmonics,Validation_Opts);
            solve_h_time = toc(solve_h_start);
        end


        %%% DEBUG
        % if iOrbit == 54
        %     debug_validation(validation_eq_terms,t0,Validation_Orbit.h,Validation_Orbit.h_dot)
        % end
        %%%

        h_analysis_start = tic;

        if Validation_Opts.get_stability
            [orbit_stab,orbit_evals] = get_h_stability(validation_eq_terms,t0);
            Validation_Orbit.evals = orbit_evals;
            % orbit_stab = h_infinite_determinant(validation_eq_terms,t0,omega,num_harmonics);
            % orbit_stab = get_h_coco_stability(Validation_Orbit,validation_eq_terms,t0,omega,num_harmonics);
        else
            orbit_stab = 1;
        end


               
        Displacement = struct("r",r,"h",Validation_Orbit.h);
        Velocity = struct("r_dot",r_dot,"h_dot",Validation_Orbit.h_dot);
        
        Eom_Terms = struct("h_inertia",h_inertia, ...
                           "h_convection",h_conv, ...
                           "h_stiffness",h_stiff, ...
                           "h_force",h_force);


        Validation_Sols(iJob) = Validation_Sols(iJob).analyse_h_solution(Displacement,Velocity,Eom_Terms,orbit_stab,Validation_Analysis_Inputs,job_orbit);
        h_analysis_time = toc(h_analysis_start);
        orbit_time = read_data_time + set_up_h_time + solve_h_time + h_analysis_time;
        job_time(iJob) = job_time(iJob) + orbit_time;

        if orbit_time < time_range(1)
            time_range(1) = orbit_time;
        end
        if orbit_time > time_range(2)
            time_range(2) = orbit_time;
        end
        fprintf("%i / %i. data: %.3f, setup: %.3f, solution: %.3f, analysis: %.3f. N_h = %i \n",...
            job_orbit,num_periodic_orbits,read_data_time,set_up_h_time,solve_h_time,h_analysis_time,num_harmonics);

        if Validation_Opts.save_orbit
            validation_name = solution_name + "\sol" + orbit_labels(job_orbit) + "_v.mat";
            save(validation_name,"-fromstruct",Validation_Orbit)
        end
    end
    time_ranges(:,iJob) = time_range;
end
mean_time = job_time./num_job_orbits;
fprintf("Mean time: %.3f, min time: %.3f, max time: %.3f \n",...
    mean(mean_time),min(time_ranges(1,:)),max(time_ranges(2,:)));

fprintf("Total time: %.3f \n",max(job_time))

Validation_Sol = combine_jobs(Validation_Sols);
end


function debug_validation(varargin) %unfortunately perpetually needed
[W_I,W_C,W_S,w_f] = varargin{1}{:};
[t0,h,h_dot] = varargin{2:end};

num_h_modes = size(varargin{1}{1},1);
num_time_points = size(t0,2);

h_ddot = zeros(num_h_modes,num_time_points);




% for iMode = 1:num_h_modes
%     h_coeff = h_frequency(iMode,:);
%     h_i = zeros(1,num_time_points) + h_coeff(1);
%     h_dot_i = zeros(1,num_time_points);
%     h_ddot_i = zeros(1,num_time_points);
%     for iHarmonic = 1:num_harmonics
%         cos_index = 1+iHarmonic;
%         sin_index = cos_index + num_harmonics;
%         cos_omega_t = cos(iHarmonic*omega*t0);
%         sin_omega_t = sin(iHarmonic*omega*t0);
% 
%         h_i = h_i + h_coeff(cos_index)*cos_omega_t + h_coeff(sin_index)*sin_omega_t;
% 
%         h_dot_i = h_dot_i + -omega*iHarmonic*h_coeff(cos_index)*sin_omega_t + omega*iHarmonic*h_coeff(sin_index)*cos_omega_t;
% 
%         h_ddot_i = h_ddot_i + -omega^2*iHarmonic^2*h_coeff(cos_index)*cos_omega_t - omega^2*iHarmonic^2*h_coeff(sin_index)*sin_omega_t;
%     end
%     h(iMode,:) = h_i;
%     h_dot(iMode,:) = h_dot_i;
%     h_ddot(iMode,:) = h_ddot_i;
% end

% h(t) and derivatives plot
% figure
% tiledlayout(num_h_modes,2)
% 
% for iMode = 1:num_h_modes
%     nexttile
%     plot(t0,h(iMode,:))
%     nexttile
%     plot(t0,h_dot(iMode,:))
%     % nexttile
%     % plot(t0,h_ddot(iMode,:))
% end



% EoM error plot
acc_force = zeros(num_h_modes,num_time_points);
conv_force = zeros(num_h_modes,num_time_points);
stiff_force = zeros(num_h_modes,num_time_points);
time_error = zeros(num_h_modes,num_time_points);
for iT = 1:num_time_points
    % acc_force(:,iT) = W_I(:,:,iT)*h_ddot(:,iT);
    conv_force(:,iT) = W_C(:,:,iT)*h_dot(:,iT);
    stiff_force(:,iT) = W_S(:,:,iT)*h(:,iT);

    % time_error(:,iT) = acc_force(:,iT) + conv_force(:,iT) + stiff_force(:,iT) -w_f(:,iT);
end

figure
tiledlayout(num_h_modes,1)
for iMode = 1:num_h_modes
    nexttile
    hold on
    % plot(t0,acc_force(iMode,:),"-r")
    plot(t0,conv_force(iMode,:),"-g")
    plot(t0,stiff_force(iMode,:),"-b")
    plot(t0,w_f(iMode,:),"y")
    % plot(t0,time_error(iMode,:),"k")
    hold off
end



%Validation equation coefficient plot
% for iQuantity = 1:4
%     h_Q = varargin{1}{iQuantity};
%     figure
%     if iQuantity < 4
%         tiledlayout(num_h_modes,num_h_modes)
%     else
%         tiledlayout(num_h_modes,1)
%     end
%     for iRow = 1:num_h_modes
%         if iQuantity < 4
%             for iCol = 1:num_h_modes
%                 nexttile
%                 plot(t0,squeeze(h_Q(iRow,iCol,:)))
%             end
%         else
%             nexttile
%             plot(t0,squeeze(h_Q(iRow,:)))
%         end
%     end
% end



end