function orbit_stab = get_h_coco_stability(Validation_Orbit,h_terms,t0,omega,num_harmonics)
[h_inertia,h_conv,h_stiff,h_force] = h_terms{:};
num_h_modes = size(h_force,1);
state_size = 2*num_h_modes;
num_time_points = size(h_force,2);

z0 = [Validation_Orbit.h;Validation_Orbit.h_dot];
% Construct A
disp_span = 1:num_h_modes;
vel_span = disp_span + num_h_modes;
system_matrix = zeros(state_size,state_size,num_time_points);
force_inertia_prod = zeros(state_size,num_time_points);
for iTime = 1:num_time_points
    system_matrix(disp_span,vel_span,iTime) = eye(num_h_modes); 
    system_matrix(vel_span,disp_span,iTime) = -h_inertia(:,:,iTime)\h_conv(:,:,iTime);
    system_matrix(vel_span,vel_span,iTime) = -h_inertia(:,:,iTime)\h_stiff(:,:,iTime);

    force_inertia_prod(vel_span,iTime) = h_inertia(:,:,iTime)\h_force(:,iTime);
end

system_matrix_frequency = time_to_frequency(system_matrix,t0,num_harmonics);
system_matrix_dt_frequency = differentiate_coefficients(system_matrix_frequency,omega);

force_frequency = time_to_frequency(force_inertia_prod,t0,num_harmonics);
force_dt_frequency = differentiate_coefficients(force_frequency,omega);

% funcs = {@(t,x,p) validation_eqation(t,x,system_matrix_frequency,force_frequency,omega)};
funcs = {@(t,x,p) validation_eqation(t,x,system_matrix_frequency,force_frequency,omega),...
            @(t,x,p) validation_eqation_dx(t,x,system_matrix_frequency,omega),...
            @(t,x,p) validation_eqation_dp(t,x),...
            @(t,x,p) validation_eqation_dt(t,x,system_matrix_dt_frequency,force_dt_frequency,omega)};

prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
%collocation Settings
coll_args = [funcs, {t0',z0', {'p'}, 0}];

prob = ode_isol2po(prob, '', coll_args{:});
bd = coco(prob, "test", [], 1, 'p', [0,0]);
eigs  = coco_bd_col(bd, {"eigs"});
orbit_stab = max(eigs);
end
%-------------------------------------------------------------------------%
%-- COCO encoding
%-------------------------------------------------------------------------%
function x_dot = validation_eqation(t,x,system_matrix_frequency,force_frequency,omega)
num_x = size(x,2);
state_size = size(x,1);


x_dot = zeros(state_size,num_x);
for iX = 1:num_x
    A = evaluate_fourier_series(system_matrix_frequency,omega,t(iX));
    force = evaluate_fourier_series(force_frequency,omega,t(iX));
    x_dot(:,iX) = A*x(:,iX) + force;
end
end
%-------------------------------------------------------------------------%
function x_dot_dx = validation_eqation_dx(t,x,system_matrix_frequency,omega)
num_x = size(x,2);
state_size = size(x,1);


x_dot_dx = zeros(state_size,state_size,num_x);
for iX = 1:num_x
    A = evaluate_fourier_series(system_matrix_frequency,omega,t(iX));
    x_dot_dx(:,:,iX) = A;
end
end
%-------------------------------------------------------------------------%
function x_dot_dp = validation_eqation_dp(t,x,p)
x_dot_dp = zeros(size(x,1),1,size(x,2));
end
%-------------------------------------------------------------------------%
function x_dot_dt = validation_eqation_dt(t,x,system_matrix_dt_frequency,force_dt_frequency,omega)
x_dot_dt = validation_eqation(t,x,system_matrix_dt_frequency,force_dt_frequency,omega);
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%-- Fourier manipulation
%-------------------------------------------------------------------------%
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
%-------------------------------------------------------------------------%
%%$ DEBUGGING
function x = evaluate_fourier_series(coeffs,omega,t)
num_dims = length(size(coeffs));
num_coefficients = size(coeffs,num_dims);
num_modes = size(coeffs,1);
coeff_dims = cell(1,num_dims);
for iDim = 1:(num_dims-1)
    coeff_dims{iDim} = 1:num_modes;
end
coeff_dims{end} = 1:num_coefficients;
num_harmonics = (num_coefficients-1)/2;


constant_index = 1;
coeff_dims{end} = constant_index;
constant_coeffs = coeffs(coeff_dims{:});

cos_index = 2:2:num_coefficients;
coeff_dims{end} = cos_index;
cos_coeffs = coeffs(coeff_dims{:});

sin_index = 3:2:num_coefficients;
coeff_dims{end} = sin_index;
sin_coeffs = coeffs(coeff_dims{:});

num_time_points = size(t,2);
x_size = size(coeffs);
x_size(end) = num_time_points;

coeff_dims{end} = 1;
x = zeros(x_size) + constant_coeffs(coeff_dims{:});
time_dims = coeff_dims;
for iHarmonic = 1:num_harmonics
    coeff_dims{end} = iHarmonic;
    for iTime = 1:num_time_points
        time_dims{end} = iTime;
        x(time_dims{:}) = x(time_dims{:}) + cos_coeffs(coeff_dims{:}).*cos(iHarmonic*omega*t(iTime));
        x(time_dims{:}) = x(time_dims{:}) + sin_coeffs(coeff_dims{:}).*sin(iHarmonic*omega*t(iTime));
    end
end
end
%------------------------------------------------------------------------%
function coeffs_dt = differentiate_coefficients(coeffs,omega)
ismatrix = ndims(coeffs) == 3;

if ismatrix
    coeffs_size = size(coeffs);
    old_size = coeffs_size([1,2]);
    coeffs = reshape(coeffs,[prod(old_size),coeffs_size(3)]);
end


num_coefficients = size(coeffs,2);
num_harmonics = (num_coefficients-1)/2;
num_modes = size(coeffs,1);

cos_index = 2:2:num_coefficients;
sin_index = 3:2:num_coefficients;

diff_prod = 1:num_harmonics;

cos_coeffs = coeffs(:,cos_index);
cos_dt_coeffs = -omega*cos_coeffs.*diff_prod;

sin_coeffs = coeffs(:,sin_index);
sin_dt_coeffs = omega*sin_coeffs.*diff_prod;

coeffs_dt = zeros(num_modes,num_coefficients);
coeffs_dt(:,cos_index) = sin_dt_coeffs;
coeffs_dt(:,sin_index) = cos_dt_coeffs;

if ismatrix
    coeffs_dt = reshape(coeffs_dt,coeffs_size);
end
end