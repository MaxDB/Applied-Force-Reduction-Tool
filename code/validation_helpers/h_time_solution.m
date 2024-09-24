function Validation_Sol = h_time_solution(Validation_Sol,Solution,Validation_Rom,solution_num)
GET_STABILITY = 0;
SAVE_ORBIT = 1;

INITIAL_NUM_HARMONICS = 5;
MAX_HARMONIC = 50;
MINIMUM_H_FORCE = 1e-6;
MAX_CONVERGENCE_ERROR = 1e-4;


num_r_modes = length(Validation_Rom.Model.reduced_modes);
disp_span = 1:num_r_modes;
vel_span = disp_span + num_r_modes;

num_h_modes = length(Validation_Rom.Dynamic_Validation_Data.current_L_modes) + num_r_modes;
h_disp_span = 1:num_h_modes;
h_vel_span = h_disp_span + num_h_modes;

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






orbit_labels = Solution.orbit_labels;
frequency = Solution.frequency;
num_periodic_orbits = length(orbit_labels);
I_L = eye(num_h_modes);

num_harmonics = INITIAL_NUM_HARMONICS;
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


    solution_converged = 0;
    while ~solution_converged
        solve_h_start = tic;
        num_time_points = length(t0);
        num_coefficients = 2*num_harmonics+1;
        B = zeros(num_time_points*num_h_modes,num_coefficients*num_h_modes);


        harmonic_span = @(index) ((index-1)*num_h_modes+1):(index*num_h_modes);
        for iTime = 1:num_time_points
            A_j = zeros(num_h_modes,num_coefficients*num_h_modes);
            A_dot_j = zeros(num_h_modes,num_coefficients*num_h_modes);
            A_ddot_j = zeros(num_h_modes,num_coefficients*num_h_modes);

            %constant terms
            A_j(:,harmonic_span(1)) = I_L;


            for nHarmonic = 1:num_harmonics
                c_nj = I_L*cos(nHarmonic*omega*t0(iTime));
                s_nj = I_L*sin(nHarmonic*omega*t0(iTime));

                term_1 = c_nj;
                term_2 = s_nj;
                A_j(:,harmonic_span(1+nHarmonic)) = term_1;
                A_j(:,harmonic_span(1+num_harmonics+nHarmonic)) = term_2;

                term_1 = -(nHarmonic*omega)*s_nj;
                term_2 = (nHarmonic*omega)*c_nj;
                A_dot_j(:,harmonic_span(1+nHarmonic)) = term_1;
                A_dot_j(:,harmonic_span(1+num_harmonics+nHarmonic)) = term_2;

                term_1 = -(nHarmonic*omega)^2*c_nj;
                term_2 = -(nHarmonic*omega)^2*s_nj;
                A_ddot_j(:,harmonic_span(1+nHarmonic)) = term_1;
                A_ddot_j(:,harmonic_span(1+num_harmonics+nHarmonic)) = term_2;
            end

            B_j = h_inertia(:,:,iTime)*A_ddot_j + h_conv(:,:,iTime)*A_dot_j + h_stiff(:,:,iTime)*A_j;

            B(harmonic_span(iTime),:) = B_j;
        end
        h_frequency_coefficients = lsqminnorm(B,reshape(h_force,num_time_points*num_h_modes,1));

        h_frequency = zeros(num_h_modes,num_coefficients);
        for iCoeff = 1:num_coefficients
            h_frequency(:,iCoeff) = h_frequency_coefficients(harmonic_span(iCoeff));
        end

        

        %-------- Check convergence -------%
        validated_h_terms = max(abs(h_force),[],2) >= MINIMUM_H_FORCE;
        if ~any(validated_h_terms)
            solution_converged = 1;
            h = frequency_to_time(h_frequency,t0,omega,num_harmonics);
            iHarmonic = num_harmonics;
        else

            h_linear = frequency_to_time(h_frequency(:,[1,2,num_harmonics+2]),t0,omega,1);
            h_n = h_linear;
            h_n_plus_two = h_n;
            harmonic_counter = 0;
            for iHarmonic = 2:num_harmonics
                harmonic_counter = harmonic_counter + 1;

                h_cos = h_frequency(:,iHarmonic+1).*cos(iHarmonic*omega*t0);
                h_sin = h_frequency(:,iHarmonic+1 + num_harmonics).*sin(iHarmonic*omega*t0);
                h_n_plus_two = h_n_plus_two + h_cos + h_sin;
                
                convergence_error = abs(h_n_plus_two - h_n)./max(abs(h_n_plus_two),[],2);
                max_error = max(convergence_error(validated_h_terms,:),[],"all");
                if harmonic_counter == 2
                    harmonic_counter = 0;
                    if max_error < MAX_CONVERGENCE_ERROR
                        solution_converged = 1;
                        
                        h = h_n_plus_two;
                        break
                    else
                        h_n = h_n_plus_two;
                        if num_harmonics > MAX_HARMONIC
                            warning("Number of harmonics exceeded maximum")
                            solution_converged = 1;
                            break
                        end
                    end
                end
            end
    
        end
        converged_harmonics = max(iHarmonic,INITIAL_NUM_HARMONICS);
        
        

        if solution_converged
            converged_harmonic_span = 2:(converged_harmonics+1);
            converged_h_frequency = h_frequency(:,[1,converged_harmonic_span,converged_harmonic_span+num_harmonics]);
            h_dot_frequency = differentiate_frequency_coefficients(converged_h_frequency,omega);
            h_dot = frequency_to_time(h_dot_frequency,t0,omega,converged_harmonics);
        end
        num_harmonics = converged_harmonics;

        solve_h_time = toc(solve_h_start);

        if ~solution_converged

            fprintf("%i / %i. solution: %.3f. Max error = %3.2e. N_h = %i \n",...
                iOrbit,num_periodic_orbits,solve_h_time,max_error,num_harmonics);
            num_harmonics = num_harmonics + 2;
            
        end


    end
    % mass = Validation_Rom.Model.mass;
    % x_dot = Validation_Rom.expand_velocity(r,r_dot,h,h_dot);
    % ke = zeros(1,num_time_points);
    % for iPoint = 1:num_time_points
    %     ke(iPoint) = 0.5*x_dot(:,iPoint)'*mass*x_dot(:,iPoint);
    % end
    h_analysis_start = tic;

    if GET_STABILITY
        orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes,num_time_points);
        for iTime = 1:(num_time_points)
            % orbit_jacobian = zeros(2*num_h_modes,2*num_h_modes);
            orbit_jacobian(h_disp_span,h_vel_span,iTime) = eye(num_h_modes);
            orbit_jacobian(h_vel_span,h_disp_span,iTime) = -h_inertia(:,:,iTime)\h_stiff(:,:,iTime);
            orbit_jacobian(h_vel_span,h_vel_span,iTime) = -h_inertia(:,:,iTime)\h_conv(:,:,iTime);
        end
        fundamental_eq = @(t,z) orbit_jacobian_func(t,z,t0,orbit_jacobian);
        [~,fundamental_mat] = ode45(fundamental_eq,[0,t0(end)],eye(2*num_h_modes));

        monodromy_mat = reshape(fundamental_mat(end,:)',2*num_h_modes,2*num_h_modes);
        orbit_evals = eig(monodromy_mat);
        orbit_stab = max(abs(orbit_evals));
    else
        orbit_stab = 1; %#ok<*UNRCH>
    end
    
    Displacement.r = r;
    Displacement.h = h;

    Velocity.r_dot = r_dot;
    Velocity.h_dot = h_dot;
    
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

    if SAVE_ORBIT
        validation_name = solution_name + "\sol" + orbit_labels(iOrbit) + "_v.mat";
        validation_orbit.h = h;
        validation_orbit.h_dot = h_dot;
        save(validation_name,"validation_orbit")
    end
end
fprintf("Mean time: %.3f, min time: %.3f, max time: %.3f \n",...
        total_time/num_periodic_orbits,time_range(1),time_range(2));
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function x_time = frequency_to_time(x_frequency,t,omega,num_harmonics)
ismatrix = ndims(x_frequency) == 3;
if ismatrix
    old_size = size(x_frequency,[1,2]);
    x_frequency = reshape(x_frequency,[prod(old_size),size(x_frequency,3)]);
    ismatrix = 1;
end

x_time = zeros(size(x_frequency,1),length(t));
for iElement = 1:size(x_time,1)
    x_time(iElement,:) = x_frequency(iElement,1);
    for iHarmonic = 1:num_harmonics
        xCos = x_frequency(iElement,1+iHarmonic)*cos(iHarmonic*omega*t);
        xSin = x_frequency(iElement,1+iHarmonic+num_harmonics)*sin(iHarmonic*omega*t);

        x_time(iElement,:) = x_time(iElement,:) + xCos + xSin;
    end

end

if ismatrix
    x_time = reshape(x_time,[old_size,length(t)]);
end
end
%-------------------------------------------------------------------------%
function x_dt_coeffs = differentiate_frequency_coefficients(x_coeffs,omega)
x_dt_coeffs = zeros(size(x_coeffs));
num_harmonics = (size(x_coeffs,2)-1)/2;
cos_span = (1:num_harmonics) + 1;
sin_span = cos_span + num_harmonics;

harmonic_mulitplier = (1:num_harmonics)*omega;
x_dt_coeffs(:,sin_span) = -harmonic_mulitplier.*x_coeffs(:,cos_span);
x_dt_coeffs(:,cos_span) = harmonic_mulitplier.*x_coeffs(:,sin_span);
end
%-------------------------------------------------------------------------%
function dz = orbit_jacobian_func(t,z,t0,orbit_jacobian)
num_h_modes = size(orbit_jacobian,1)/2;
approx_jacobian = zeros(2*num_h_modes);
for iRow = 1:(2*num_h_modes)
    for iCol = 1:(2*num_h_modes)
        approx_jacobian(iRow,iCol) = interp1(t0,squeeze(orbit_jacobian(iRow,iCol,:)),t);
    end
end

dz = approx_jacobian*reshape(z,2*num_h_modes,2*num_h_modes);
dz = reshape(dz,4*num_h_modes^2,1);
end
