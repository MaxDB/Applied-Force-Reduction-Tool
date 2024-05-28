function Dyn_Data = h_time_solution(Dyn_Data,Validation_Rom,solution_num)
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

% Validation_Input = Validation_Rom.get_solver_inputs("h_prediction_test");
% h_terms = @(r,r_dot,r_ddot) get_h_error_terms_test(r,r_dot,r_ddot,Validation_Input);

Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("h_analysis");

solution_name = Rom.data_path + "dynamic_sol_" + solution_num;
sol_labels = Dyn_Data.solution_labels{1,solution_num};
frequency = Dyn_Data.frequency{1,solution_num};
num_periodic_orbits = length(sol_labels);
I_L = eye(num_h_modes);

Dyn_Data = Dyn_Data.pre_allocate_validation_data(solution_num,num_periodic_orbits);
num_harmonics = INITIAL_NUM_HARMONICS;
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
    

    h_analysis_start = tic;
    Dyn_Data = Dyn_Data.analyse_h_solution(r,r_dot,h,h_dot,Validation_Analysis_Inputs,solution_num,iOrbit);
    h_analysis_time = toc(h_analysis_start);

    fprintf("%i / %i. data: %.3f, setup: %.3f, solution: %.3f, analysis: %.3f. N_h = %i \n",...
        iOrbit,num_periodic_orbits,read_data_time,set_up_h_time,solve_h_time,h_analysis_time,num_harmonics);
end

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

