function [ke_mode_tilde,ke_condensed_tilde,ke_mode_hat,ke_condensed_hat] = h_kinetic_energy(r,r_dot,h,h_dot,Eom_Input)
num_x = size(r,2);

num_coupling_coeffs = size(Eom_Input.Disp_Data.g_L_coeffs,2);
num_h_coupling_coeffs = size(Eom_Input.H_Coupling_Gradient.beta_bar,1);
num_coeffs = size(Eom_Input.input_order,1);


scale_factor = Eom_Input.scale_factor;
shift_factor = Eom_Input.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r - shift_factor);

theta_H_beta_bar = Eom_Input.Disp_Data.theta_H_beta_bar;
G_beta_bar = Eom_Input.H_Coupling_Gradient.beta_bar;
G_theta_H_beta = Eom_Input.beta_G_theta_H;
theta_H_G_beta = Eom_Input.beta_theta_H_G;

num_r_modes = size(r,1);
num_h_modes = size(h,1);
num_L_modes = num_h_modes - num_r_modes;

h_r_span = 1:num_r_modes;
h_L_span = (1:num_L_modes) + num_r_modes;


ke_mode_tilde = zeros(num_h_modes,num_x);
ke_condensed_tilde = zeros(1,num_x);
ke_mode_hat = zeros(num_h_modes,num_x);
ke_condensed_hat = zeros(1,num_x);

input_order = Eom_Input.input_order;
for iX = 1:num_x
    r_transformed_i = r_transformed(:,iX);
    r_dot_i = r_dot(:,iX);

    h_i = h(:,iX);
    h_dot_i = h_dot(:,iX);

    % for iTerm = 1:num_coeffs
    %     r_power_products(iTerm,1) = prod(r_transformed_i.^input_order(:,iTerm));
    % end
    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_r_modes
        r_power_products = r_power_products.*r_transformed_i(iMode).^input_order(:,iMode);
    end

    r_products_theta = r_power_products(1:num_coupling_coeffs,:);
    r_dr_products_theta = r_products_theta(Eom_Input.Disp_Data.diff_mapping{1,1}).*Eom_Input.Disp_Data.diff_scale_factor{1,1};

    r_products_G = r_power_products(1:num_h_coupling_coeffs,:);
    r_dr_products_G = r_products_G(Eom_Input.H_Coupling_Gradient.diff_mapping{1,1}).*Eom_Input.H_Coupling_Gradient.diff_scale_factor{1,1};
    
    %%% r modes
    r_dot_i_hat = r_dot_i+h_dot_i(h_r_span);
    ke_r_tilde = 1/2*(r_dot_i.*r_dot_i);
    ke_r_hat = 1/2*(r_dot_i_hat.*r_dot_i_hat);
    
    %%% L modes
    g_dot_L = Eom_Input.Disp_Data.g_L_coeffs*r_dr_products_theta*r_dot_i;
    L_dot_i_hat = g_dot_L + h_dot_i(h_L_span);
    ke_L_tilde = 1/2*(g_dot_L.*g_dot_L);
    ke_L_hat = 1/2*(L_dot_i_hat.*L_dot_i_hat);

    %%% H modes
    theta_dr_prod = r_dr_products_theta*r_dot_i;
    G_dr_prod = (r_dr_products_G*r_dot_i);
    % theta_dr_H * theta_Dr_H
    theta_dr_theta_dr_prod = theta_dr_prod'*theta_H_beta_bar*theta_dr_prod;
    % theta_dr_H * h_coupling_dr
    pre_theta_dr_G_dr_prod = tensorprod(theta_H_G_beta,G_dr_prod,3,1);
    theta_dr_G_dr_prod = theta_dr_prod'*pre_theta_dr_G_dr_prod*h_i;
    % theta_dr_H * h_coupling
    pre_theta_dr_G_prod = tensorprod(theta_H_G_beta,r_products_G,3,1);
    theta_dr_G_prod = theta_dr_prod'*pre_theta_dr_G_prod*h_dot_i;
    % h_coupling_dr * theta_dr_H
    pre_G_dr_theta_dr_prod = squeeze(tensorprod(G_dr_prod',G_theta_H_beta,2,1));
    G_dr_theta_dr_prod = h_i'*pre_G_dr_theta_dr_prod*theta_dr_prod;
    % h_coupling_dr * h_coupling_dr
    pre_G_dr_G_dr_prod = squeeze(tensorprod(G_dr_prod',G_beta_bar,2,1));
    G_dr_G_dr_prod = h_i'*tensorprod(pre_G_dr_G_dr_prod,G_dr_prod,3,1)*h_i;
    % h_coupling_dr * h_coupling
    G_dr_G_prod = h_i'*tensorprod(pre_G_dr_G_dr_prod,r_products_G,3,1)*h_dot_i;
    % h_coulping * theta_dr_H
    pre_G_theta_dr_prod = squeeze(tensorprod(r_products_G',G_theta_H_beta,2,1));
    G_theta_dr_prod = h_dot_i'*pre_G_theta_dr_prod*theta_dr_prod;
    % h_coupling * h_coupling_dr
    pre_G_G_dr_prod = squeeze(tensorprod(r_products_G',G_beta_bar,2,1));
    G_G_dr_prod = h_dot_i'*tensorprod(pre_G_G_dr_prod,G_dr_prod,3,1)*h_i;
    % h_coupling * h_coupling
    G_G_prod = h_dot_i'*tensorprod(pre_G_G_dr_prod,r_products_G,3,1)*h_dot_i;

    
    ke_H_tilde = 1/2*theta_dr_theta_dr_prod;
    ke_H_hat = 1/2*(theta_dr_theta_dr_prod+theta_dr_G_dr_prod+G_dr_theta_dr_prod+G_dr_G_dr_prod...
        + theta_dr_G_prod + G_dr_theta_dr_prod + G_dr_G_prod + G_theta_dr_prod + G_G_dr_prod + G_G_prod);
    
    %%%%
    ke_mode_tilde(:,iX) = [ke_r_tilde;ke_L_tilde];
    ke_condensed_tilde(1,iX) = ke_H_tilde;
    ke_mode_hat(:,iX) = [ke_r_hat;ke_L_hat];
    ke_condensed_hat(1,iX) = ke_H_hat;

end
end