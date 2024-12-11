function Validation_Sol = solve_h_prediction(Validation_Sol,Solution,Validation_Rom,solution_num)
Validation_Opts = Validation_Sol.Validation_Options;

num_r_modes = length(Validation_Rom.Model.reduced_modes);
disp_span = 1:num_r_modes;
vel_span = disp_span + num_r_modes;


Solution_Type = Solution.Solution_Type;
orbit_type = Solution_Type.orbit_type;
solution_name = Validation_Rom.data_path + "dynamic_sol_" + solution_num;

total_time = 0;
time_range = [inf,0];
%%% Set up h-problem

switch orbit_type
    case "free"
        Eom_Input = Validation_Rom.get_solver_inputs("coco_backbone");
        reduced_eom = @(t,z,zeta) coco_eom(t,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
        %
        Validation_Input = Validation_Rom.get_solver_inputs("h_prediction");
        h_terms = @(r,r_dot,r_ddot) get_h_error_terms(r,r_dot,r_ddot,Validation_Input);

        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("h_analysis");
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

num_harmonics = Validation_Opts.initial_harmonic;
for iOrbit = 1:num_periodic_orbits
    read_data_start = tic;
    sol = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(iOrbit));
    read_data_time = toc(read_data_start);

    set_up_h_start = tic;
    t0 = sol.tbp';
    x = sol.xbp';
    r = x(disp_span,:);
    r_dot = x(vel_span,:);
    omega = frequency(1,iOrbit);

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
    end
    %h_inertia * h_ddot  +  h_conv * h_dot  +  h_stiff * h  =  h_force
    set_up_h_time = toc(set_up_h_start);

    validation_eq_terms = {h_inertia,h_conv,h_stiff,h_force};
    r_force = Validation_Rom.Force_Polynomial.evaluate_polynomial(r);
    solution_converged = 0;
    while ~solution_converged
        solve_h_start = tic;
        h_frequency = h_solver(validation_eq_terms,t0,omega,num_harmonics);
        [solution_converged,num_harmonics,Validation_Orbit] = check_h_convergence(validation_eq_terms,r_force,h_frequency,t0,omega,num_harmonics,Validation_Opts);
        solve_h_time = toc(solve_h_start);
    end

    h_analysis_start = tic;

    if Validation_Opts.get_stability
        orbit_stab = get_h_stability(validation_eq_terms,t0);
        % orbit_stab = h_infinite_determinant(validation_eq_terms,t0,omega,num_harmonics);
        % orbit_stab = get_h_coco_stability(Validation_Orbit,validation_eq_terms,t0,omega,num_harmonics);
    else
        orbit_stab = 1;
    end

    Displacement.r = r;
    Displacement.h = Validation_Orbit.h;

    Velocity.r_dot = r_dot;
    Velocity.h_dot = Validation_Orbit.h_dot;

    Eom_Terms.h_inertia = h_inertia;
    Eom_Terms.h_convection = h_conv;
    Eom_Terms.h_stiffness = h_stiff;
    Eom_Terms.h_force = h_force;

    Validation_Sol = Validation_Sol.analyse_h_solution(Displacement,Velocity,Eom_Terms,orbit_stab,Validation_Analysis_Inputs,iOrbit);
    h_analysis_time = toc(h_analysis_start);
    orbit_time = read_data_time + set_up_h_time + solve_h_time + h_analysis_time;
    total_time = total_time + orbit_time;

    if orbit_time < time_range(1)
        time_range(1) = orbit_time;
    end
    if orbit_time > time_range(2)
        time_range(2) = orbit_time;
    end
    fprintf("%i / %i. data: %.3f, setup: %.3f, solution: %.3f, analysis: %.3f. N_h = %i \n",...
        iOrbit,num_periodic_orbits,read_data_time,set_up_h_time,solve_h_time,h_analysis_time,num_harmonics);

    if Validation_Opts.save_orbit
        validation_name = solution_name + "\sol" + orbit_labels(iOrbit) + "_v.mat";
        save(validation_name,"Validation_Orbit")
    end
end
fprintf("Mean time: %.3f, min time: %.3f, max time: %.3f \n",...
    total_time/num_periodic_orbits,time_range(1),time_range(2));
end


