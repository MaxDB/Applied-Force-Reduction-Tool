function [ke_tilde,ke_hat] = h_kinetic_energy(r,r_dot,h,h_dot,Eom_Input)
num_x = size(r,2);
num_r_modes = size(r,1);

num_disp_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp_r_disp,3);
num_grad_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp_r_disp,1);
num_coeffs = size(Eom_Input.input_order,1);


scale_factor = Eom_Input.scale_factor;
shift_factor = Eom_Input.shift_factor;

r_transformed = scale_factor.*(r + shift_factor);

h_disp_beta_bar = Eom_Input.Beta_Bar_Data.h_disp;
h_disp_r_disp_beta_bar = Eom_Input.Beta_Bar_Data.h_disp_r_disp;
r_disp_h_disp_beta_bar = permute(h_disp_r_disp_beta_bar,[3,2,1]);
r_disp_beta_bar = Eom_Input.Beta_Bar_Data.r_disp;

input_order = Eom_Input.input_order;

ke_tilde = zeros(1,num_x);
ke_hat = zeros(1,num_x);
for iX = 1:num_x
    r_dot_i = r_dot(:,iX);
    r_transformed_i = r_transformed(:,iX);
    
    h_i = h(:,iX);
    h_dot_i = h_dot(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_r_modes
        r_power_products = r_power_products.*r_transformed_i(iMode).^input_order(:,iMode);
    end
    

    r_products_disp = r_power_products(1:num_disp_coeffs,:);
    r_dr_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,1}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,1};

    r_products_grad = r_power_products(1:num_grad_coeffs);
    r_dr_products_grad = r_products_grad(Eom_Input.Disp_Grad_Data.diff_mapping{1,1}).*Eom_Input.Disp_Grad_Data.diff_scale_factor{1,1};


    %--
    r_dr_r_dot_prod = r_dr_products_disp*r_dot_i;
    r_dr_grad_r_dot_prod = r_dr_products_grad*r_dot_i;
    
    %--
    r_disp_r_disp = r_dr_r_dot_prod'*r_disp_beta_bar*r_dr_r_dot_prod;
    
    %--
    r_disp_h_disp_prod = r_dr_r_dot_prod'*tensorprod(r_disp_h_disp_beta_bar,r_dr_grad_r_dot_prod,3,1);
    r_disp_h_disp = r_disp_h_disp_prod*h_i;

    %--
    r_disp_h_vel_prod = r_dr_r_dot_prod'*tensorprod(r_disp_h_disp_beta_bar,r_products_grad,3,1);
    r_disp_h_vel = r_disp_h_vel_prod*h_dot_i;

    %--
    h_disp_prod = tensorprod(r_dr_grad_r_dot_prod',h_disp_beta_bar,2,1);
    h_disp_prod = tensorprod(h_i',h_disp_prod,2,2);
    h_disp_h_disp_prod = tensorprod(h_disp_prod,r_dr_grad_r_dot_prod,4,1);
    h_disp_h_disp = tensorprod(h_disp_h_disp_prod,h_i,3,1);

    %--
    h_disp_h_vel_prod = tensorprod(h_disp_prod,r_products_grad,4,1);
    h_disp_h_vel = tensorprod(h_disp_h_vel_prod,h_dot_i,3,1);

    %--
    h_vel_prod = tensorprod(r_products_grad',h_disp_beta_bar,2,1);
    h_vel_prod = tensorprod(h_dot_i',h_vel_prod,2,2);
    h_vel_h_vel_prod = tensorprod(h_vel_prod,r_products_grad,4,1);
    h_vel_h_vel = tensorprod(h_vel_h_vel_prod,h_dot_i,3,1);

    %--

    ke_tilde(:,iX) = 0.5*r_disp_r_disp;
    ke_hat(:,iX) = 0.5*(r_disp_r_disp + h_disp_h_disp + h_vel_h_vel) + r_disp_h_disp + r_disp_h_vel + h_disp_h_vel;

end
end