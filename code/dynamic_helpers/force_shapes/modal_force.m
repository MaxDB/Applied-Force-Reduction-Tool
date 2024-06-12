function force = modal_force(time,amplitude,period,mode_map)
num_modes = size(mode_map,2);
num_x = size(time,2);
force_amplitude = zeros(num_modes,num_x);
force_amplitude(mode_map,:) = amplitude;
force = force_amplitude.*sinpi(time*2./period);
end