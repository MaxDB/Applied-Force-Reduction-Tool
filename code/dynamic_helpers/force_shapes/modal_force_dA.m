function force_dA = modal_force_dA(time,amplitude,period,mode_map)
num_modes = size(mode_map,2);
num_forces = size(amplitude,1);
num_x = size(time,2);
force_dA = zeros(num_modes,num_forces,num_x);
mode_index = find(mode_map);
for iForce = 1:num_forces
    force_dA(mode_index(iForce),iForce,:) = sinpi(time*2./period);
end
% force = force_amplitude.*sinpi(time*2/period);
end