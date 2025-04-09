function x_dot = direct_forced_eom(~,x,modal_restoring_force,damping)
num_x = size(x,2);
num_modes = size(x,1)/2;

disp_span = 1:num_modes;
q = x(disp_span,:);

vel_span = disp_span + num_modes;
q_dot = x(vel_span,:);

x_dot = zeros(2*num_modes,num_x);
x_dot(disp_span,:) = q_dot;
for iX = 1:num_x
    x_dot(vel_span,iX) = -modal_restoring_force(q(:,iX)) - damping*q_dot(:,iX);
end
end