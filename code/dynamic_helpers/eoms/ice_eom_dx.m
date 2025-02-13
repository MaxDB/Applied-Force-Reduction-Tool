function x_dot_dx = ice_eom_dx(~,x,zeta,input_order,Force_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

num_force_coeffs = size(Force_Data.coeffs,2);
num_coeffs = size(input_order,1);

disp_span = 1:num_modes;
r = x(disp_span,:);
scale_factor = Force_Data.scale_factor;
shift_factor = Force_Data.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r - shift_factor);

vel_span = disp_span + num_modes;

x_dot_dx = zeros(2*num_modes,2*num_modes,num_x);


I_R = eye(num_modes);
for iX = 1:num_x
    r_i = r_transformed(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
    end
    
    r_products_force = r_power_products(1:num_force_coeffs,:);
    r_dr_products_force = r_products_force(Force_Data.diff_mapping{1,1}).*Force_Data.diff_scale_factor{1,1};

    %-------------
    restoring_force_dr = Force_Data.coeffs*r_dr_products_force;
    
    %-------------
    x_dot_dx(disp_span,vel_span,iX) = I_R;
    x_dot_dx(vel_span,disp_span,iX) = -restoring_force_dr;
    x_dot_dx(vel_span,vel_span,iX) = -zeta(iX)*I_R;
end

end