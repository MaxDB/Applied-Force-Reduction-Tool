function force = sine_force(time,amplitude,period)
force = amplitude.*sinpi(time*2./period);
end