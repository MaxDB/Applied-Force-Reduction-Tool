function force_dA = sine_force_dA(time,amplitude,period)
amp_size = size(amplitude,1);
force_dA = eye(amp_size).*sinpi(time*2/period);
end