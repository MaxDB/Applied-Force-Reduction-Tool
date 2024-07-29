function x_dot = coco_eom(~,x,zeta,input_order,Force_Data,Disp_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

num_force_coeffs = size(Force_Data.coeffs,2);
num_coupling_coeffs = size(Disp_Data.beta_bar,1);
num_coeffs = size(input_order,1);

disp_span = 1:num_modes;
r = x(disp_span,:);
scale_factor = Force_Data.scale_factor;
shift_factor = Force_Data.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r - shift_factor);


vel_span = disp_span + num_modes;
r_dot = x(vel_span,:);

x_dot = zeros(2*num_modes,num_x);
x_dot(disp_span,:) = r_dot;

% r_power_products = ones(num_coeffs,1);


for iX = 1:num_x
    r_i = r_transformed(:,iX);
    r_dot_i = r_dot(:,iX);
    % for iTerm = 1:num_coeffs
    %     r_power_products(iTerm,1) = prod(r_transformed(:,iX).^input_order(:,iTerm));
    % end
    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
    end

    r_products_coupling = r_power_products(1:num_coupling_coeffs,:);
    restoring_force = Force_Data.coeffs*r_power_products(1:num_force_coeffs,:);

    r_dr_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,2}).*Disp_Data.diff_scale_factor{1,2};
    
    

    theta_prod = r_dr_products_coupling'*Disp_Data.beta_bar;
    
    inertia_term = eye(num_modes) + theta_prod*r_dr_products_coupling;
    convection_term = theta_prod*(tensorprod(r_dr2_products_coupling,r_dot_i,3,1)*r_dot_i);

    x_dot(vel_span,iX) = -inertia_term\(convection_term+restoring_force) - zeta(iX).*r_dot_i;
end
end