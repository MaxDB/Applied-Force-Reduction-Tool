function orbit_stability = h_infinite_determinant(h_terms,t0,omega,num_harmonics)
[h_inertia,h_conv,h_stiff,h_force] = h_terms{:};
num_h_modes = size(h_force,1);
state_size = 2*num_h_modes;
num_time_points = size(h_force,2);



% % TEST
% 
% num_harmonics = 1;
% state_size = 2;
% omega = 2;
% 
% delta = 0.4; %1
% epsilon = 0.4;
% 
% A0 = [0,1;-delta,0];
% B = [0,0;-2*epsilon,0];
% 
% t0 = linspace(0,2*pi/omega,num_time_points);
% system_matrix = zeros(state_size,state_size,num_time_points);
% for iTime = 1:num_time_points
%     system_matrix(:,:,iTime) = A0 + B*cos(omega*t0(iTime));
% end

% Construct A
disp_span = 1:num_h_modes;
vel_span = disp_span + num_h_modes;
system_matrix = zeros(state_size,state_size,num_time_points);
for iTime = 1:num_time_points
    system_matrix(disp_span,vel_span,iTime) = eye(num_h_modes); 
    system_matrix(vel_span,disp_span,iTime) = -h_inertia(:,:,iTime)\h_conv(:,:,iTime);
    system_matrix(vel_span,vel_span,iTime) = -h_inertia(:,:,iTime)\h_stiff(:,:,iTime);
end
%%%%%%%%%%%%%%%%%%%%
% %normalise inputs
% T = 2*pi/omega;
% t0 = t0/T;
% omega = 2*pi;



% 
% scale_matrix = zeros(num_h_modes);
% for iMode = 1:num_h_modes
%     diagonal_component = squeeze(system_matrix(vel_span(iMode),vel_span(iMode),:));
%     diagonal_amplitude = std(diagonal_component);
%     scale_matrix(iMode,iMode) = 1/diagonal_amplitude;
% end
% 
% system_matrix(vel_span,:,:)  = pagemtimes(scale_matrix,system_matrix(vel_span,:,:));
%%%%%%%%%%%%%%%%%%%

% harmonic components of A
system_matrix_frequency = time_to_frequency(system_matrix,t0,num_harmonics);
%wasted time on zero and indentiy component



% construct Hill's matrix
num_coefficients = 2*num_harmonics + 1;
hill_matrix = zeros(num_coefficients*state_size,num_coefficients*state_size);
diagonal_counter = 0;
for iHarmonic = -num_harmonics:num_harmonics
    for iState = 1:state_size
        diagonal_counter = diagonal_counter + 1;
        hill_matrix(diagonal_counter,diagonal_counter) = -1j*iHarmonic*omega;
    end
end

row_counter = 0;
col_range = @(index) (index:state_size:(num_coefficients*state_size));
for k_harmonic = -num_harmonics:num_harmonics
   
    
    for iMode = 1:state_size %h mode
        row_counter = row_counter + 1;
        for jMode = 1:state_size %% stiffness columns
            C_m = c_coeffs(squeeze(system_matrix_frequency(iMode,jMode,:)),k_harmonic,num_harmonics);
            col_indices = col_range(jMode);
            hill_matrix(row_counter,col_indices) = hill_matrix(row_counter,col_indices) + C_m.';
        end
    end
    
end

truncated_harmonic = 3;


characteristic_exponents = eig(hill_matrix);
unstable_eig = characteristic_exponents(real(characteristic_exponents) >= 0);
characteristic_multipliers = exp(characteristic_exponents*(2*pi)/omega);

orbit_stability = max(abs(characteristic_multipliers));

%%% TEST
hill_dim = num_coefficients*state_size;
orbit_stability = 1+max(real(characteristic_exponents));
TOL = 1e-3;
for iEig = size(unstable_eig,1)

    alg_multiplicity = nnz(ismembertol(imag(unstable_eig),imag(unstable_eig(iEig)),TOL));
    geo_multiplicity = hill_dim - rank(hill_matrix - unstable_eig(iEig)*eye(hill_dim),TOL);

    if alg_multiplicity > geo_multiplicity
        orbit_stability = 1000;
        break
    end
end
end

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
harmonic_span = -num_harmonics:num_harmonics;
for iElement = 1:num_elements
    X = fourier_coefficients(t0,x_time(iElement,:));
    X_sym = fftshift(X);
    centre_point = floor(size(X_sym,2)/2) + 1;
    coeffs = X_sym(centre_point+harmonic_span);

    x_frequency(iElement,:) = coeffs;
end

if ismatrix
    x_frequency = reshape(x_frequency,[old_size,num_coeffs]);
end
end
%-------------------------------------------------------------------------%
function X = fourier_coefficients(t,x)
% T = t(end)-t(1);   %sample interval

% dtMin = T/length(t);
% tLin = t(1):dtMin:t(end);
max_points = length(t);
num_points = ceil(max_points*1);
t_lin = linspace(t(1),t(end),num_points);

L = length(t_lin);   %discrete length

x_lin = interp1(t,x,t_lin);

% x_signal = [x_lin(1),repmat(x_lin(2:end),1,5)];
x_signal = x_lin(1:(end-1)); %!!!!!!

X = fft(x_signal);
X = X/(L-1);
end
%-------------------------------------------------------------------------%
%-- Harmonic Balance Helpers
%-------------------------------------------------------------------------%
function C_m = c_coeffs(d,k_harmonic,num_harmonics)
num_coefficients = 2*num_harmonics + 1;
centre_point = floor(num_coefficients/2) + 1;

C_m = zeros(num_coefficients,1);

upper_limit = min(num_harmonics,num_harmonics + k_harmonic);
lower_limit = max(-num_harmonics,-num_harmonics + k_harmonic);

for m_harmonic = lower_limit:upper_limit
    C_m(m_harmonic + centre_point) = C_m(m_harmonic + centre_point) + d(k_harmonic-m_harmonic + centre_point);
end
end
%-------------------------------------------------------------------------%
