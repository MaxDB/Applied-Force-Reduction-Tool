function x_dot_dt = direct_forced_eom_dt(t,x,amp,period,~,~,modal_applied_force)
num_x = size(x,2);
num_modes = size(x,1)/2;

vel_span = (1:num_modes) + num_modes;

x_dot_dt = zeros(2*num_modes,num_x);
x_dot_dt(vel_span,:) = amp*modal_applied_force.*(2*pi./period).*cos(2*pi./period.*t);
end