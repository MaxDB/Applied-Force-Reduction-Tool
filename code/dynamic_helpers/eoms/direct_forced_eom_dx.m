function x_dot_dx = direct_forced_eom_dx(~,x,~,~,modal_stiffness,modal_damping,~)
num_x = size(x,2);
num_modes = size(x,1)/2;

disp_span = 1:num_modes;
q = x(disp_span,:);

vel_span = disp_span + num_modes;

x_dot_dx = zeros(2*num_modes,2*num_modes,num_x);
I_N = eye(num_modes);

for iX = 1:num_x
    r_cell = num2cell(q(:,iX));
    dfdx = modal_stiffness(r_cell{:});

    x_dot_dx(disp_span,vel_span,iX) = I_N;
    x_dot_dx(vel_span,vel_span,iX) = -modal_damping;
    x_dot_dx(vel_span,disp_span,iX) = -dfdx;
end

end