function h_frequency_alt = h_harmonic_balance(h_terms,t0,omega,num_harmonics)

[h_inertia,h_conv,h_stiff,h_force] = h_terms{:};
num_h_modes = size(h_force,1);

h_force_frequency = time_to_frequency(h_force,t0,num_harmonics);
h_inertia_frequency = time_to_frequency(h_inertia,t0,num_harmonics);
h_conv_frequency = time_to_frequency(h_conv,t0,num_harmonics);
h_stiff_frequency = time_to_frequency(h_stiff,t0,num_harmonics);


%------------------------- TEST ----------------------%
% max_points = length(t0);
% num_points = ceil(max_points*1);
% t_lin = linspace(t0(1),t0(end),num_points);
% 
% 
% if num_harmonics == 23
%     test_name = "h_stiff";
%     test = eval(test_name);
%     test_frequency = eval(test_name + "_frequency");
%     test_time = evaluate_fourier_series(test_frequency,omega,t0);
%     num_h_modes = size(test_time,1);
%     figure;
%     tiledlayout("flow")
%     for iPlot = 1:num_h_modes
%         for jPlot = 1:num_h_modes
%             nexttile
%             hold on
%             plot(t0,squeeze(test(iPlot,jPlot,:)))
%             plot(t0,squeeze(test_time(iPlot,jPlot,:)),"--")
% 
%             x_lin = interp1(t0,squeeze(test(iPlot,jPlot,:)),t_lin);
%             plot(t_lin,x_lin,"x")
%             hold off
%         end
%     end
%     a
% end
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

%convert to [cos, sin] ceofficient form
h_frequency_alt = zeros(size(h_frequency));
alt_cos_index = (1:num_harmonics) + 1;
alt_sin_index = alt_cos_index + num_harmonics;

h_frequency_alt(:,constant_index) = h_frequency(:,constant_index);
h_frequency_alt(:,alt_cos_index) = h_frequency(:,cos_index);
h_frequency_alt(:,alt_sin_index) = h_frequency(:,sin_index);
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
function X = fourier_coefficients(t,x)
% T = t(end)-t(1);   %sample interval

% dtMin = T/length(t);
% tLin = t(1):dtMin:t(end);
max_points = length(t);
num_points = ceil(max_points*0.8);
t_lin = linspace(t(1),t(end),num_points);

L = length(t_lin);   %discrete length

x_lin = interp1(t,x,t_lin);

% x_signal = [x_lin(1),repmat(x_lin(2:end),1,5)];
x_signal = x_lin(1:(end-1));

X = fft(x_signal);
X = X/(L-1);
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