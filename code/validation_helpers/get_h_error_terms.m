function  [h_inertia,h_conv,h_stiff,h_force] = get_h_error_terms(r,r_dot,r_ddot,Eom_Input)
num_x = size(r,2);
num_r_modes = size(r,1);

num_reduced_force_coeffs = size(Eom_Input.Reduced_Force_Data.coeffs,2);
num_coupling_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp_r_disp,3);
num_h_force_grad_coeffs = size(Eom_Input.H_Force_Data.coeffs,3);
num_h_coupling_grad_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp,1);
num_coeffs = size(Eom_Input.input_order,1);

num_h_modes = size(Eom_Input.H_Force_Data.coeffs,2);

scale_factor = Eom_Input.scale_factor;
shift_factor = Eom_Input.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r + shift_factor);

h_disp_beta_bar = Eom_Input.Beta_Bar_Data.h_disp;
h_disp_r_disp_beta_bar = Eom_Input.Beta_Bar_Data.h_disp_r_disp;

input_order = Eom_Input.input_order;


h_inertia = zeros(num_h_modes,num_h_modes,num_x);
h_conv = zeros(num_h_modes,num_h_modes,num_x);
h_stiff = zeros(num_h_modes,num_h_modes,num_x);
h_force = zeros(num_h_modes,num_x);
for iX = 1:num_x
    r_ddot_i = r_ddot(:,iX);
    r_dot_i = r_dot(:,iX);
    r_transformed_i = r_transformed(:,iX);


    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_r_modes
        r_power_products = r_power_products.*r_transformed_i(iMode).^input_order(:,iMode);
    end
    
    %---
    r_force = Eom_Input.Reduced_Force_Data.coeffs*r_power_products(1:num_reduced_force_coeffs,:);
    h_force_grad = tensorprod(Eom_Input.H_Force_Data.coeffs,r_power_products(1:num_h_force_grad_coeffs,:),3,1);
    
    %--
    r_products_disp = r_power_products(1:num_coupling_coeffs,:);
    r_dr_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,1}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,2}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,2};


    %%% Inertia
    H_beta_prod = squeeze(tensorprod(r_products_disp',h_disp_beta_bar,2,1));
    h_inertia(:,:,iX) = tensorprod(H_beta_prod,r_products_disp,3,1);

    %%% Convection
    r_dr_r_dot_prod = r_dr_products_disp*r_dot_i;
    H_r_dr_r_dot_prod = tensorprod(H_beta_prod,r_dr_r_dot_prod,3,1);
    h_conv(:,:,iX) = 2*H_r_dr_r_dot_prod;

    %%% Stiffness
    r_d2r_r_dot_r_dot_prod = tensorprod(r_dr2_products_disp,r_dot_i,3,1)*r_dot_i;
    r_dr_r_ddot_prod = r_dr_products_disp*r_ddot_i;
    stiff_sum = r_d2r_r_dot_r_dot_prod + r_dr_r_ddot_prod;
    stiff_prod = tensorprod(H_beta_prod,stiff_sum,3,1);

    h_stiff(:,:,iX) = stiff_prod + h_force_grad;

    %%% Force
    H_beta_disp_prod = squeeze(tensorprod(r_products_disp',h_disp_r_disp_beta_bar,2,1));
    force_1 = H_beta_disp_prod*stiff_sum;
    force_1(1:num_r_modes,1) = force_1(1:num_r_modes,1) + r_force;
    h_force(:,iX) = -force_1;
end
end

