function force_dT = sine_force_dTper(time,amplitude,period)
force_dT = amplitude.*cospi(2*time./period).*(-2*pi*time./period.^2);
end