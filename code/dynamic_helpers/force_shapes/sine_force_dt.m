function force_dt = sine_force_dt(time,amplitude,period)
force_dt = amplitude.*cospi(time*2./period).*2*pi./period;
end