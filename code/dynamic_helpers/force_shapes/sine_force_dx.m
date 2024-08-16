function force_dx = sine_force_dx(time,amplitude,period,num_r_modes)
amp_size = size(amplitude,1);
num_points = size(time,2);
force_dx = zeros(amp_size,num_r_modes*2,num_points);
end