function force_dT = modal_force_dTper(time,amplitude,period,mode_map)
num_modes = size(mode_map,2);
num_x = size(time,2);
force_dT = zeros(num_modes,num_x);
force_dT(mode_map,:) = amplitude.*cospi(2*time./period).*(-2*pi*time./period.^2);
% force = force_amplitude.*sinpi(time*2/period);
end