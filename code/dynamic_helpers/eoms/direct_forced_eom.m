function x_dot = direct_forced_eom(t,x,amp,period,modal_restoring_force,modal_damping,modal_applied_force)
num_x = size(x,2);
num_modes = size(x,1)/2;

disp_span = 1:num_modes;
q = x(disp_span,:);

vel_span = disp_span + num_modes;
q_dot = x(vel_span,:);

x_dot = zeros(2*num_modes,num_x);
x_dot(disp_span,:) = q_dot;

frequency = 2*pi./period;
sin_arg = frequency.*t;
for iX = 1:num_x
    q_i = q(:,iX);
    q_dot_i = q_dot(:,iX);
    sin_arg_i = sin_arg(iX);

    %--
    restoring_force = modal_restoring_force(q_i);
    damping_force = modal_damping*q_dot_i;
    applied_force = amp*modal_applied_force*sin(sin_arg_i);

    %--
    x_dot(vel_span,iX) = applied_force - damping_force - restoring_force;
end
end