function  [h_inertia,h_conv,h_stiff,h_force] = get_h_error_terms(r,r_dot,r_ddot,Eom_Input)
num_x = size(r,2);
num_r_modes = size(r,1);

num_reduced_force_coeffs = size(Eom_Input.Reduced_Force_Data.coeffs,2);
num_disp_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp_r_disp,3);
num_h_force_grad_coeffs = size(Eom_Input.H_Force_Data.coeffs,3);
num_h_disp_grad_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp,1);
num_coeffs = size(Eom_Input.input_order,1);

num_h_modes = size(Eom_Input.H_Force_Data.coeffs,2);

scale_factor = Eom_Input.scale_factor;
shift_factor = Eom_Input.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r + shift_factor);

h_disp_beta_bar = Eom_Input.Beta_Bar_Data.h_disp;
h_disp_r_disp_beta_bar = Eom_Input.Beta_Bar_Data.h_disp_r_disp;

input_order = Eom_Input.input_order;

if num_h_disp_grad_coeffs > num_disp_coeffs
    disp_diff_mapping = Eom_Input.Disp_Grad_Data.diff_mapping;
    disp_diff_scale_factor = Eom_Input.Disp_Grad_Data.diff_scale_factor;
    num_max_disp_coeffs = num_h_disp_grad_coeffs;
else
    disp_diff_mapping = Eom_Input.Physical_Disp_Data.diff_mapping;
    disp_diff_scale_factor = Eom_Input.Physical_Disp_Data.diff_scale_factor;
    num_max_disp_coeffs = num_disp_coeffs;
end


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
    % r_products_disp_all = r_power_products(1:num_max_disp_coeffs,:);
    % 
    % dr_dr_power_products = r_products_disp_all(disp_diff_mapping{1,1}).*disp_diff_scale_factor{1,1};
    % d2r_dr2_products =  r_products_disp_all(disp_diff_mapping{1,2}).*disp_diff_scale_factor{1,2};
    % 
    % 
    % r_dr_products_disp = dr_dr_power_products(1:num_disp_coeffs,:);
    % r_dr2_products_disp = d2r_dr2_products(1:num_disp_coeffs,:,:);
    % 
    % r_products_grad = r_products_disp_all(1:num_h_disp_grad_coeffs,:);
    % r_dr_products_grad = dr_dr_power_products(1:num_h_disp_grad_coeffs,:);
    % r_dr2_products_grad = d2r_dr2_products(1:num_h_disp_grad_coeffs,:,:);

    r_products_disp = r_power_products(1:num_disp_coeffs,:);
    r_dr_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,1}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,2}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,2};


    r_products_grad = r_power_products(1:num_h_disp_grad_coeffs);
    r_dr_products_grad = r_products_grad(Eom_Input.Disp_Grad_Data.diff_mapping{1,1}).*Eom_Input.Disp_Grad_Data.diff_scale_factor{1,1};
    r_dr2_products_grad = r_products_grad(Eom_Input.Disp_Grad_Data.diff_mapping{1,2}).*Eom_Input.Disp_Grad_Data.diff_scale_factor{1,2};

    % r_products_grad = r_products_disp;
    % r_dr_products_grad = r_dr_products_disp;
    % r_dr2_products_grad = r_dr2_products_disp;

    %%% Inertia
    H_beta_prod = squeeze(tensorprod(r_products_grad',h_disp_beta_bar,2,1));
    h_inertia(:,:,iX) = tensorprod(H_beta_prod,r_products_grad,3,1);

    %%% Convection
    r_dr_r_dot_prod = r_dr_products_grad*r_dot_i;
    H_r_dr_r_dot_prod = tensorprod(H_beta_prod,r_dr_r_dot_prod,3,1);
    h_conv(:,:,iX) = 2*H_r_dr_r_dot_prod;

    %%% Stiffness
    r_d2r_r_dot_r_dot_prod = tensorprod(r_dr2_products_grad,r_dot_i,3,1)*r_dot_i;
    r_dr_r_ddot_prod = r_dr_products_grad*r_ddot_i;
    stiff_sum = r_d2r_r_dot_r_dot_prod + r_dr_r_ddot_prod;
    stiff_prod = tensorprod(H_beta_prod,stiff_sum,3,1);

    h_stiff(:,:,iX) = stiff_prod + h_force_grad;

    %%% Force
    H_beta_disp_prod = squeeze(tensorprod(r_products_grad',h_disp_r_disp_beta_bar,2,1));
    disp_r_d2r_r_dot_r_dot_prod = tensorprod(r_dr2_products_disp,r_dot_i,3,1)*r_dot_i;
    disp_r_dr_r_ddot_prod = r_dr_products_disp*r_ddot_i;
    disp_stiff_sum = disp_r_d2r_r_dot_r_dot_prod + disp_r_dr_r_ddot_prod;
    force_1 = H_beta_disp_prod*disp_stiff_sum;
    force_1(1:num_r_modes,1) = force_1(1:num_r_modes,1) + r_force;

    % force_1(1:num_r_modes,1) = 0;
    h_force(:,iX) = -force_1;
end
end

