function x_dot_dzeta = direct_eom_dzeta(~,x)
num_x = size(x,2);
num_modes = size(x,1)/2;

vel_span = (1:num_modes) + num_modes;
q_dot = x(vel_span,:);


x_dot_dzeta = zeros(2*num_modes,1,num_x);
x_dot_dzeta(vel_span,1,:) = -q_dot;
end