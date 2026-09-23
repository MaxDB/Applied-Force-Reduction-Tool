function x_dot_dTper = direct_forced_eom_dTper(t,x,amp,period,~,~,modal_applied_force,Harmonic_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

vel_span = (1:num_modes) + num_modes;

x_dot_dTper = zeros(2*num_modes,1,num_x);

if isempty(Harmonic_Data)
    force_time_dTper = (-2*pi*t./period.^2).*cos(2*pi./period.*t);
else
    force_time_dTper = Harmonic_Data.harmonics_dT(t,2*pi./period);
end

x_dot_dTper(vel_span,1,:) = amp*modal_applied_force.*force_time_dTper;
end