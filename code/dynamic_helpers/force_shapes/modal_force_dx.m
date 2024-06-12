function force_dx = modal_force_dx(time,amplitude,period,mode_map)
num_modes = size(mode_map,2);
num_x = size(time,2);
force_dx = zeros(num_modes,num_modes*2,num_x);
end