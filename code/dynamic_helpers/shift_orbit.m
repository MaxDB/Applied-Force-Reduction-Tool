function [t_shift,x_shift] = shift_orbit(t,x)
    period = t(end) - t(1);
    omega = 2*pi/period;

    x_frequency = time_to_frequency(x,t,1);
    a1 = x_frequency(2);
    b1 = x_frequency(3);
    %find phase diff compared to cosine

    phase = atan2(-b1,a1);
    t_shifted = t + phase/omega;
    t_all = [t_shifted(1:(end-1))-period,t_shifted,t_shifted(2:end)+period];
    x_all = [x(:,(1:(end-1))),x,x(:,2:end)];

    keep_indices = t_all >= 0 & t_all <=period; 
    keep_indices(find(keep_indices,1)-1) = 1;
    keep_indices(find(keep_indices,1,"last")+1) = 1;
    
    t_shift = t_all(keep_indices);
    x_shift = x_all(keep_indices);
end

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
num_points = ceil(max_points*0.9);
t_lin = linspace(t(1),t(end),num_points);

% L = length(tLin);   %discrete length

x_lin = interp1(t,x,t_lin);
X = fft(x_lin);
X = X/num_points;
end
