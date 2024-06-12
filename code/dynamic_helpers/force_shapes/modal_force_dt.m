function force_dt = modal_force_dt(time,amplitude,period,mode_map)
num_modes = size(mode_map,2);
num_x = size(time,2);
force_amplitude = zeros(num_modes,num_x);
force_amplitude(mode_map,:) = amplitude;
force_dt = force_amplitude.*cospi(time*2./period).*2*pi./period;
end