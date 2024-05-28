function Dyn_Data = h_harmonic_balance(Dyn_Data,Validation_Rom,solution_num)
INITIAL_NUM_HARMONICS = 5;
MINIMUM_H_FORCE = 1e-6;
MAX_CONVERGENCE_ERROR = 1e-4;

num_r_modes = length(Validation_Rom.Model.reduced_modes);
disp_span = 1:num_r_modes;
vel_span = disp_span + num_r_modes;

num_h_modes = length(Validation_Rom.Dynamic_Validation_Data.current_L_modes) + num_r_modes;

Rom = Dyn_Data.Dynamic_Model;
%%% Set up h-problem
Eom_Input = Rom.get_solver_inputs("coco_backbone");
reduced_eom = @(t,z,zeta) coco_eom(t,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
% 
Validation_Input = Validation_Rom.get_solver_inputs("h_prediction");
h_terms = @(r,r_dot,r_ddot) get_h_error_terms(r,r_dot,r_ddot,Validation_Input);

% Validation_Input_Test = Validation_Rom.get_solver_inputs("h_prediction_test");
% h_terms = @(r,r_dot,r_ddot) get_h_error_terms_test(r,r_dot,r_ddot,Validation_Input_Test);

solution_name = Rom.data_path + "dynamic_sol_" + solution_num;
sol_labels = Dyn_Data.solution_labels{1,solution_num};
frequency = Dyn_Data.frequency{1,solution_num};
num_periodic_orbits = length(sol_labels);

num_harmonics = INITIAL_NUM_HARMONICS;
Dyn_Data = Dyn_Data.pre_allocate_validation_data(solution_num,num_periodic_orbits);
for iOrbit = 1:num_periodic_orbits
    read_data_start = tic;
    sol = po_read_solution('',convertStringsToChars(solution_name),sol_labels(iOrbit));
    read_data_time = toc(read_data_start);
    
    set_up_h_start = tic;
    t0 = sol.tbp';
    x = sol.xbp';
    r = x(disp_span,:);
    r_dot = x(vel_span,:);
    omega = frequency(1,iOrbit);

    x_dot = reduced_eom(t0,x,zeros(size(t0)));
    r_ddot  = x_dot(vel_span,:);

    [h_inertia,h_conv,h_stiff,h_force] = h_terms(r,r_dot,r_ddot); %lots of scope to speed up
    %h_inertia * h_ddot  +  h_conv * h_dot  +  h_stiff * h  =  h_force
    set_up_h_time = toc(set_up_h_start);


    
    solution_converged = 0;
    while ~solution_converged
        solve_hb_start = tic;
        h_force_frequency = time_to_frequency(h_force,t0,num_harmonics);
        h_inertia_frequency = time_to_frequency(h_inertia,t0,num_harmonics);
        h_conv_frequency = time_to_frequency(h_conv,t0,num_harmonics);
        h_stiff_frequency = time_to_frequency(h_stiff,t0,num_harmonics);


        %------------------------- TEST ----------------------%
        % P_test = frequency_to_time(P_frequency,t0,omega);
        % inertia_test = frequency_to_time(h_inertia_frequency,t0,omega);
        % stiffness_test = frequency_to_time(h_stiffness_frequency,t0,omega);

        % test_harmonics(t0,P,P_test)
        % test_harmonics(t0,h_inertia,inertia_test)
        % test_harmonics(t0,h_stiffness,stiffness_test)

        % num_harmonics = 1;
        % t_test = linspace(0,t0(end),200);
        % x_time = sin(omega*t_test);
        % x_frequency = time_to_frequency(x_time,t_test);
        %
        % x_test = frequency_to_time(x_frequency,t_test,omega);
        %
        % test_harmonics(t_test,x_time,x_test)
        %------------------------- TEST ----------------------%
        num_coefficients = size(h_force_frequency,2);
        force_frequency_coeffs = h_force_frequency';

        c_stiffness = permute(h_stiff_frequency,[3,1,2]);
        c_conv = permute(h_conv_frequency,[3,1,2]);
        c_inertia = permute(h_inertia_frequency,[3,1,2]);

        constant_index = 1;
        cos_index = 2:2:num_coefficients;
        sin_index = 3:2:num_coefficients;

        ai_0 = c_inertia(constant_index,:,:);
        ai_n = c_inertia(cos_index,:,:);
        bi_n = c_inertia(sin_index,:,:);

        ac_0 = c_conv(constant_index,:,:);
        ac_n = c_conv(cos_index,:,:);
        bc_n = c_conv(sin_index,:,:);

        as_0 = c_stiffness(constant_index,:,:);
        as_n = c_stiffness(cos_index,:,:);
        bs_n = c_stiffness(sin_index,:,:);

        upsilon = zeros(num_h_modes*num_coefficients,num_h_modes*num_coefficients);

        % go through equation (D.1) working out the coefficients for A and B for each harmonic

        % constant term
        row_counter = 0;
        col_range = @(index) ((index-1)*num_coefficients + 1):(index*num_coefficients);
        for iMode = 1:num_h_modes %h mode

            % c_0
            row_counter = row_counter + 1;
            for jMode = 1:num_h_modes %% stiffness columns
                C_n = c_0(as_0(1,iMode,jMode),as_n(:,iMode,jMode),num_harmonics);
                col_indices = col_range(jMode);
                upsilon(row_counter,col_indices) = C_n';
            end


            for kHarmonic = 1:num_harmonics
                %cos terms
                row_counter = row_counter + 1;
                %inertia terms
                inertia_sf = - (kHarmonic*omega)^2;
                for jMode = 1:num_h_modes
                    C_n = c_k(ai_0(1,iMode,jMode),ai_n(:,iMode,jMode),bi_n(:,iMode,jMode),kHarmonic,num_harmonics);
                    C_n(1) = 0; %no constant term
                    col_indices = col_range(jMode);
                    upsilon(row_counter,col_indices) = upsilon(row_counter,col_indices) + inertia_sf*C_n';
                end

                %convective terms
                convective_sf = (kHarmonic*omega);
                for jMode = 1:num_h_modes
                    C_n = d_k(ac_0(1,iMode,jMode),ac_n(:,iMode,jMode),bc_n(:,iMode,jMode),kHarmonic,num_harmonics);
                    C_n(1) = 0; %no constant term
                    col_indices = col_range(jMode);
                    upsilon(row_counter,col_indices) = upsilon(row_counter,col_indices) + convective_sf*C_n';
                end

                %stiffness terms
                % c_k
                for jMode = 1:num_h_modes
                    C_n = c_k(as_0(1,iMode,jMode),as_n(:,iMode,jMode),bs_n(:,iMode,jMode),kHarmonic,num_harmonics);
                    col_indices = col_range(jMode);
                    upsilon(row_counter,col_indices) = upsilon(row_counter,col_indices) + C_n';
                end

                % sin terms
                row_counter = row_counter + 1;
                %inertia terms
                inertia_sf = - (kHarmonic*omega)^2;
                for jMode = 1:num_h_modes
                    C_n = d_k(ai_0(1,iMode,jMode),ai_n(:,iMode,jMode),bi_n(:,iMode,jMode),kHarmonic,num_harmonics);
                    C_n(1) = 0; %no constant term
                    col_indices = col_range(jMode);
                    upsilon(row_counter,col_indices) = upsilon(row_counter,col_indices) + inertia_sf*C_n';
                end

                %convective terms
                convective_sf = -(kHarmonic*omega);
                for jMode = 1:num_h_modes
                    C_n = c_k(ac_0(1,iMode,jMode),ac_n(:,iMode,jMode),bc_n(:,iMode,jMode),kHarmonic,num_harmonics);
                    C_n(1) = 0; %no constant term
                    col_indices = col_range(jMode);
                    upsilon(row_counter,col_indices) = upsilon(row_counter,col_indices) + convective_sf*C_n';
                end


                %stiffness terms
                % d_k
                for jMode = 1:num_h_modes
                    C_n = d_k(as_0(1,iMode,jMode),as_n(:,iMode,jMode),bs_n(:,iMode,jMode),kHarmonic,num_harmonics);
                    col_indices = col_range(jMode);
                    upsilon(row_counter,col_indices) = upsilon(row_counter,col_indices) + C_n';
                end
            end

        end



        force_frequency_coeffs = reshape(force_frequency_coeffs,[num_h_modes*num_coefficients,1]);
        h_frequency_coeffs = upsilon\force_frequency_coeffs;

        h_frequency = reshape(h_frequency_coeffs,num_coefficients,num_h_modes)';
        % h = frequency_to_time(h_frequency,t0,omega,num_harmonics);
        

        %-------- Check convergence -------%
        validated_h_terms = max(abs(h_force),[],2) > MINIMUM_H_FORCE;
        if ~any(validated_h_terms)
            solution_converged = 1;
            h = frequency_to_time(h_frequency,t0,omega,num_harmonics);
        else

            h_linear = frequency_to_time(h_frequency(:,1:3),t0,omega,1);
            h_n = h_linear;
            h_n_plus_two = h_n;
            harmonic_counter = 0;
            for iHarmonic = 2:num_harmonics
                harmonic_counter = harmonic_counter + 1;

                h_cos = h_frequency(:,iHarmonic*2).*cos(iHarmonic*omega*t0);
                h_sin = h_frequency(:,iHarmonic*2+1).*sin(iHarmonic*omega*t0);
                h_n_plus_two = h_n_plus_two + h_cos + h_sin;
                convergence_error = abs(h_n_plus_two - h_n)./max(abs(h_n_plus_two),[],2);
                max_error = max(convergence_error(validated_h_terms,:),[],"all");
                if harmonic_counter == 2
                    harmonic_counter = 0;
                    if max_error < MAX_CONVERGENCE_ERROR
                        solution_converged = 1;
                        num_harmonics = max(iHarmonic,INITIAL_NUM_HARMONICS);
                        h = h_n_plus_two;
                        break
                    else
                        h_n = h_n_plus_two;
                    end
                end
            end

        end
        
        solve_hb_time = toc(solve_hb_start);

        if ~solution_converged

            fprintf("%i / %i. solution: %.3f. Max error = %3.2e. N_h = %i \n",...
                iOrbit,num_periodic_orbits,solve_hb_time,max_error,num_harmonics);
            num_harmonics = num_harmonics + 2;
        end

    end

    h_analysis_start = tic;
    Dyn_Data = Dyn_Data.analyse_h_solution(r,h,Validation_Rom,solution_num,iOrbit);
    h_analysis_time = toc(h_analysis_start);

    fprintf("%i / %i. data: %.3f, setup: %.3f, solution: %.3f, analysis: %.3f. N_h = %i \n",...
        iOrbit,num_periodic_orbits,read_data_time,set_up_h_time,solve_hb_time,h_analysis_time,num_harmonics);
end
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%-- Fourier manipulation
%-------------------------------------------------------------------------%
function x_frequency = time_to_frequency(x_time,t0,num_harmonics)


ismatrix = ndims(x_time) == 3;

if ismatrix
    x_size = size(x_time);
    old_size = x_size([1,2]);
    x_time = reshape(x_time,[prod(old_size),x_size(3)]);
end

num_coeffs = 2*num_harmonics+1;
num_elements = size(x_time,1);
x_frequency = zeros(num_elements,num_coeffs);
for iElement = 1:num_elements
    X = fourier_coefficients(t0,x_time(iElement,:));
    alpha0 = real(X(1,1));
    alpha = 2*real(X(1,2:(num_harmonics+1)));
    beta = -2*imag(X(1,2:(num_harmonics+1)));

    coeffs = zeros(2*num_harmonics+1,1);
    coeffs(1) = alpha0;
    for iHarmonic = 1:num_harmonics
        coeffs(iHarmonic*2) = alpha(iHarmonic);
        coeffs(iHarmonic*2+1) = beta(iHarmonic);
    end
    x_frequency(iElement,:) = coeffs;
end

if ismatrix
    x_frequency = reshape(x_frequency,[old_size,num_coeffs]);
end
end
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
        xCos = x_frequency(iElement,iHarmonic*2)*cos(iHarmonic*omega*t);
        xSin = x_frequency(iElement,iHarmonic*2+1)*sin(iHarmonic*omega*t);

        x_time(iElement,:) = x_time(iElement,:) + xCos + xSin;
    end

end

if ismatrix
    x_time = reshape(x_time,[old_size,length(t)]);
end
end
%-------------------------------------------------------------------------%
function X = fourier_coefficients(t,x)
% T = t(end)-t(1);   %sample interval

% dtMin = T/length(t);
% tLin = t(1):dtMin:t(end);
max_points = length(t);
num_points = ceil(max_points*0.9);
t_lin = linspace(t(1),t(end),num_points);

% L = length(tLin);   %discrete length

x_lin = interp1(t,x,t_lin);
X = fft(x_lin);
X = X/num_points;
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%-- Harmonic Balance Helpers
%-------------------------------------------------------------------------%
function C_n = c_0(a_0,a_n,num_harmonics)
A_n = zeros(num_harmonics,1);
B_n = zeros(num_harmonics,1);

A_0 = a_0(1);
for iHarmonic = 1:num_harmonics
    A_n(iHarmonic) = A_n(iHarmonic) + 0.5*a_n(iHarmonic);
    B_n(iHarmonic) = B_n(iHarmonic) + 0.5*a_n(iHarmonic); % dont "correct"
end
C_n = convert_to_C(A_0,A_n,B_n,num_harmonics);
end
%-------------------------------------------------------------------------%
function C_n = c_k(a_0,a_n,b_n,kHarmonic,num_harmonics)
A_n = zeros(num_harmonics,1);
B_n = zeros(num_harmonics,1);

A_n(kHarmonic) = A_n(kHarmonic) + a_0;
A_0 = a_n(kHarmonic);

for nHarmonic = 1:(kHarmonic-1)
    A_n(nHarmonic) = A_n(nHarmonic) + 0.5*a_n(kHarmonic-nHarmonic);
    B_n(nHarmonic) = B_n(nHarmonic) - 0.5*b_n(kHarmonic-nHarmonic);
end

for nHarmonic = 1:(num_harmonics-kHarmonic)
    A_n(nHarmonic) = A_n(nHarmonic) + 0.5*a_n(nHarmonic+kHarmonic);
    B_n(nHarmonic) = B_n(nHarmonic) + 0.5*b_n(nHarmonic+kHarmonic);
end

for nHarmonic = (kHarmonic+1):(num_harmonics)
    A_n(nHarmonic) = A_n(nHarmonic) + 0.5*a_n(nHarmonic-kHarmonic);
    B_n(nHarmonic) = B_n(nHarmonic) + 0.5*b_n(nHarmonic-kHarmonic);
end

C_n = convert_to_C(A_0,A_n,B_n,num_harmonics);
end
%-------------------------------------------------------------------------%
function C_n = d_k(a_0,a_n,b_n,kHarmonic,num_harmonics)
A_n = zeros(num_harmonics,1);
B_n = zeros(num_harmonics,1);

B_n(kHarmonic) = B_n(kHarmonic) + a_0;
A_0 = B_n(kHarmonic);

for nHarmonic = 1:(kHarmonic-1)
    A_n(nHarmonic) = A_n(nHarmonic) + 0.5*b_n(kHarmonic-nHarmonic);
    B_n(nHarmonic) = B_n(nHarmonic) + 0.5*a_n(kHarmonic-nHarmonic);
end

for nHarmonic = 1:(num_harmonics-kHarmonic)
    A_n(nHarmonic) = A_n(nHarmonic) + 0.5*b_n(nHarmonic+kHarmonic);
    B_n(nHarmonic) = B_n(nHarmonic) - 0.5*a_n(nHarmonic+kHarmonic);
end

for nHarmonic = (kHarmonic+1):(num_harmonics)
    A_n(nHarmonic) = A_n(nHarmonic) - 0.5*b_n(nHarmonic-kHarmonic);
    B_n(nHarmonic) = B_n(nHarmonic) + 0.5*a_n(nHarmonic-kHarmonic);
end

C_n = convert_to_C(A_0,A_n,B_n,num_harmonics);
end
%-------------------------------------------------------------------------%
function C_n = convert_to_C(A_0,A_n,B_n,num_harmonics)
num_coeffs = 2*num_harmonics+1;
C_n = zeros(num_coeffs,1);

cos_index = 2:2:num_coeffs;
sin_index = 3:2:num_coeffs;
C_n(1) = A_0;
C_n(cos_index) = A_n;
C_n(sin_index) = B_n;
end
%-------------------------------------------------------------------------%