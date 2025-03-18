function x_dot = ice_eom(~,x,zeta,input_order,Force_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

num_force_coeffs = size(Force_Data.coeffs,2);
num_coeffs = size(input_order,1);

disp_span = 1:num_modes;
r = x(disp_span,:);
scale_factor = Force_Data.scale_factor;
shift_factor = Force_Data.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r + shift_factor);


vel_span = disp_span + num_modes;
r_dot = x(vel_span,:);

x_dot = zeros(2*num_modes,num_x);
x_dot(disp_span,:) = r_dot;
 


for iX = 1:num_x
    r_i = r_transformed(:,iX);
    r_dot_i = r_dot(:,iX);
    
    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
    end


    restoring_force = Force_Data.coeffs*r_power_products(1:num_force_coeffs,:);
    x_dot(vel_span,iX) = -restoring_force - zeta(iX).*r_dot_i;
end
end