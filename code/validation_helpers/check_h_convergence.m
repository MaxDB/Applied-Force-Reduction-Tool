function [solution_converged,num_harmonics,Validation_Orbit] = check_h_convergence(h_terms,r_force,h_frequency,t0,omega,num_harmonics,Validation_Opts)
[h_inertia,h_conv,h_stiff,h_force] = h_terms{:};
max_r_force = max(abs(r_force),[],2);
num_r_modes = size(max_r_force,1);
solution_converged = 0;
validated_h_terms = max(abs(h_force/max(max_r_force)),[],2) >= Validation_Opts.minimum_validation_force;
% validated_h_terms(1:num_r_modes) = true;
if ~any(validated_h_terms)
    solution_converged = 1;
    h = frequency_to_time(h_frequency,t0,omega,num_harmonics);
    iHarmonic = num_harmonics;
else

    h_dot_frequency = differentiate_frequency_coefficients(h_frequency,omega);
    h_ddot_frequency = differentiate_frequency_coefficients(h_dot_frequency,omega);
    
    h_linear = frequency_to_time(h_frequency(:,[1,2,num_harmonics+2]),t0,omega,1);
    h_dot_linear = frequency_to_time(h_dot_frequency(:,[1,2,num_harmonics+2]),t0,omega,1);
    h_ddot_linear = frequency_to_time(h_ddot_frequency(:,[1,2,num_harmonics+2]),t0,omega,1);



    h_n = h_linear;
    h_dot_n = h_dot_linear;
    h_ddot_n = h_ddot_linear;

    h_n_plus_two = h_n;
    h_dot_n_plus_two = h_dot_n;
    h_ddot_n_plus_two = h_ddot_n;
    harmonic_counter = 0;
    
    eom_error_n = get_eom_error(h_n,h_dot_n,h_ddot_n,h_inertia,h_conv,h_stiff,h_force,0);
    

    for iHarmonic = 2:num_harmonics
        harmonic_counter = harmonic_counter + 1;
        
        cos_t = cos(iHarmonic*omega*t0);
        sin_t = sin(iHarmonic*omega*t0);

        h_cos = h_frequency(:,iHarmonic+1).*cos_t;
        h_sin = h_frequency(:,iHarmonic+1 + num_harmonics).*sin_t;
        h_n_plus_two = h_n_plus_two + h_cos + h_sin;

        h_dot_sin = h_dot_frequency(:,iHarmonic+1).*sin_t;
        h_dot_cos = h_dot_frequency(:,iHarmonic+1 + num_harmonics).*cos_t;
        h_dot_n_plus_two = h_dot_n_plus_two + h_dot_cos + h_dot_sin;

        h_ddot_cos = h_ddot_frequency(:,iHarmonic+1).*cos_t;
        h_ddot_sin = h_ddot_frequency(:,iHarmonic+1 + num_harmonics).*sin_t;
        h_ddot_n_plus_two = h_ddot_n_plus_two + h_ddot_cos + h_ddot_sin;

        eom_error_n_plus_two = get_eom_error(h_n_plus_two,h_dot_n_plus_two,h_ddot_n_plus_two,h_inertia,h_conv,h_stiff,h_force,0);
        
        convergence_error = get_convergence_error(eom_error_n,eom_error_n_plus_two);
        % convergence_error = get_convergence_error(h_n,h_n_plus_two);
        convergence_error = eom_error_n_plus_two;


        [max_row_error,max_error_row_index] = max(convergence_error(validated_h_terms,:),[],2);
        [max_error,max_error_index] = max(max_row_error);
        
        % if num_harmonics > 50
        %     hold on
        %     plot(t0,h_n(12,:),"--")
        %     plot(t0,h_n_plus_two(12,:),"-")
        %     hold off
        % end

        if harmonic_counter == 2
            harmonic_counter = 0;
            if max_error < Validation_Opts.maximum_convergence_error
                solution_converged = 1;

                h = h_n_plus_two;
          
                break
            else
                if iHarmonic > Validation_Opts.maximum_harmonic
                    warning("Number of harmonics exceeded maximum")
                    solution_converged = 1;
                    h = h_n;
 
                    break
                end
                 h_n = h_n_plus_two;
                 h_dot_n = h_dot_n_plus_two;
                 h_ddot_n = h_ddot_n_plus_two;
                 eom_error_n = eom_error_n_plus_two;
            end
        end
    end

end
converged_harmonics = max(iHarmonic,Validation_Opts.initial_harmonic);



if solution_converged
    converged_harmonic_span = 2:(converged_harmonics+1);
    converged_h_frequency = h_frequency(:,[1,converged_harmonic_span,converged_harmonic_span+num_harmonics]);
    h_dot_frequency = differentiate_frequency_coefficients(converged_h_frequency,omega);
    h_dot = frequency_to_time(h_dot_frequency,t0,omega,converged_harmonics);
    
    
    Validation_Orbit.h = h;
    Validation_Orbit.h_dot = h_dot;


    % ---------------------------------%
    % DEBUG
    % h_ddot_frequency = differentiate_frequency_coefficients(h_dot_frequency,omega);
    % h_ddot = frequency_to_time(h_ddot_frequency,t0,omega,converged_harmonics);
    % 
    % num_h_modes = size(h,1);
    % num_time_points = size(t0,2);
    % inertia_force = zeros(num_h_modes,num_time_points);
    % convective_force = zeros(num_h_modes,num_time_points);
    % restoring_force = zeros(num_h_modes,num_time_points);
    % for iTime = 1:num_time_points
    %     inertia_force(:,iTime) = h_inertia(:,:,iTime)*h_ddot(:,iTime);
    %     convective_force(:,iTime) = h_conv(:,:,iTime)*h_dot(:,iTime);
    %     restoring_force(:,iTime) = h_stiff(:,:,iTime)*h(:,iTime);
    % end
    % net_force = inertia_force + convective_force + restoring_force - h_force;
    % figure;
    % tiledlayout("flow")
    % for iMode = 1:num_h_modes
    %     nexttile
    %     hold on
    %     plot(t0,inertia_force(iMode,:),"b")
    %     plot(t0,convective_force(iMode,:),"r")
    %     plot(t0,restoring_force(iMode,:),"g")
    %     plot(t0,h_force(iMode,:),"y")
    %     plot(t0,net_force(iMode,:),"k")
    %     hold off
    % end
    % ---------------------------------%
end
num_harmonics = converged_harmonics;



if ~solution_converged
    fprintf("Max error = %3.2e. N_h = %i \n",...
        max_error,num_harmonics);
    num_harmonics = num_harmonics + 2;
    Validation_Orbit = [];
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
%---------------------
function convergence_error = get_convergence_error(x_1,x_2)
convergence_error = abs(x_2 - x_1)./max(abs(x_2),[],2);
end
%--------------------
function eom_error = get_eom_error(h,h_dot,h_ddot,h_inertia,h_conv,h_stiff,h_force,debug)
num_t_points = size(h,2);
num_h_modes = size(h,1);

eom_error = zeros(num_h_modes,num_t_points);
norm_eom_error = zeros(num_h_modes,num_t_points);
acc_force = zeros(num_h_modes,num_t_points);
stiff_force = zeros(num_h_modes,num_t_points);
conv_force = zeros(num_h_modes,num_t_points);
for iT = 1:num_t_points
    acc_force(:,iT) = h_inertia(:,:,iT)*h_ddot(:,iT);
    conv_force(:,iT) = h_conv(:,:,iT)*h_dot(:,iT);
    stiff_force(:,iT) = h_stiff(:,:,iT)*h(:,iT);

    max_force = max(abs([acc_force(:,iT),conv_force(:,iT),stiff_force(:,iT),h_force(:,iT)]),[],2);

    eom_error(:,iT) = (acc_force(:,iT) + conv_force(:,iT) + stiff_force(:,iT) - h_force(:,iT));
    norm_eom_error(:,iT) = eom_error(:,iT)./max_force;
end

if ~debug
    return
end
figure
tiledlayout(num_h_modes,1);
for iMode = 1:num_h_modes
    nexttile
    hold on
    plot(acc_force(iMode,:),"r")
    plot(conv_force(iMode,:),"b")
    plot(stiff_force(iMode,:),"g")
    plot(h_force(iMode,:),"m")
    plot(eom_error(iMode,:),"k")
    hold off
end
end
%-------------------
