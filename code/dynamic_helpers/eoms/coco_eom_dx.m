function x_dot_dx = coco_eom_dx(~,x,zeta,input_order,Force_Data,Disp_Data)
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

r_transformed = scale_factor.*(r + shift_factor);

vel_span = disp_span + num_modes;
r_dot = x(vel_span,:);

x_dot_dx = zeros(2*num_modes,2*num_modes,num_x);

switch num_modes
    case 1
        d2_dims = 2;
    otherwise
        d2_dims = 3;
end


I_R = eye(num_modes);
for iX = 1:num_x
    r_i = r_transformed(:,iX);
    r_dot_i = r_dot(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
    end
    
    r_products_force = r_power_products(1:num_force_coeffs,:);
    r_dr_products_force = r_products_force(Force_Data.diff_mapping{1,1}).*Force_Data.diff_scale_factor{1,1};
    
    r_products_coupling = r_power_products(1:num_coupling_coeffs,:);
    r_dr_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,2}).*Disp_Data.diff_scale_factor{1,2};
    r_dr3_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,3}).*Disp_Data.diff_scale_factor{1,3};
    
    %-------------
    disp_dr_prod = r_dr_products_coupling'*Disp_Data.beta_bar;
    disp_dr2_prod = tensorprod(pagetranspose(r_dr2_products_coupling),Disp_Data.beta_bar,2,1);

    inertia = disp_dr_prod*r_dr_products_coupling;
    inertia_dr = tensorprod(disp_dr_prod,r_dr2_products_coupling,2,1) ...
        + tensorprod(disp_dr2_prod,r_dr_products_coupling,d2_dims,1);
    
    %-------------
    r_dr2_r_dot_prod = tensorprod(r_dr2_products_coupling,r_dot_i,3,1);
    r_dr3_r_dot_prod = tensorprod(r_dr3_products_coupling,r_dot_i,4,1);
    
    convection_dr_dot = disp_dr_prod*(r_dr2_r_dot_prod);
    pre_convection_dr = tensorprod(disp_dr2_prod,r_dr2_r_dot_prod,d2_dims,1) ...
        + tensorprod(disp_dr_prod,r_dr3_r_dot_prod,2,1);
    convection_dr = tensorprod(pre_convection_dr,r_dot_i,3,1);
    convection = convection_dr_dot*r_dot_i;

    %-------------
    restoring_force = Force_Data.coeffs*r_products_force;
    restoring_force_dr = Force_Data.coeffs*r_dr_products_force;
    
    %-------------
    r_ddot = -inertia\(convection+restoring_force);
    pre_r_ddot_dr = tensorprod(inertia_dr,r_ddot,3,1);
    pre_r_ddot_dr_dot = (convection_dr + restoring_force_dr);
    %-------------
    x_dot_dx(disp_span,vel_span,iX) = I_R;
    x_dot_dx(vel_span,disp_span,iX) = -inertia\(pre_r_ddot_dr+pre_r_ddot_dr_dot);
    x_dot_dx(vel_span,vel_span,iX) = -inertia\(2*convection_dr_dot) - zeta(iX)*I_R;
end

end