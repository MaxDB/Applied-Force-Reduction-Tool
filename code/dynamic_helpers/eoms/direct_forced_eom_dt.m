function x_dot_dt = direct_forced_eom_dt(t,x,amp,period,~,~,modal_applied_force,Harmonic_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

vel_span = (1:num_modes) + num_modes;

x_dot_dt = zeros(2*num_modes,num_x);

if isempty(Harmonic_Data)
    force_time_dt = (2*pi./period).*cos(2*pi./period.*t);
else
    force_time_dt = Harmonic_Data.harmonics_dt(t,2*pi./period);
end

x_dot_dt(vel_span,:) = amp*modal_applied_force.*force_time_dt;
end