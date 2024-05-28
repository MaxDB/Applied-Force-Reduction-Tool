function  [h_inertia,h_conv,h_stiff,h_force] = get_h_error_terms(r,r_dot,r_ddot,Eom_Input)
num_x = size(r,2);

num_reduced_force_coeffs = size(Eom_Input.Reduced_Force_Data.coeffs,2);
num_coupling_coeffs = size(Eom_Input.Condensed_Disp_Data.beta_L,2);
num_h_stiffness_coeffs = size(Eom_Input.H_Stiffness_Data.coeffs,3);
num_h_coupling_coeffs = size(Eom_Input.H_Coupling_Gradient.beta_bar,1);
num_coeffs = size(Eom_Input.input_order,1);


scale_factor = Eom_Input.scale_factor;
shift_factor = Eom_Input.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r - shift_factor);

G_beta_bar = Eom_Input.H_Coupling_Gradient.beta_bar;
L_theta_beta = Eom_Input.Condensed_Disp_Data.beta_L;
G_theta_beta = Eom_Input.beta_G_theta;



input_order = Eom_Input.input_order';
num_h_modes = size(Eom_Input.H_Stiffness_Data.coeffs,1);
I_h = eye(num_h_modes);
r_power_products = ones(num_coeffs,1);

switch num_h_modes
    case 1
        d2_dims = 2;
    otherwise
        d2_dims = 3;
end

h_inertia = zeros(num_h_modes,num_h_modes,num_x);
h_conv = zeros(num_h_modes,num_h_modes,num_x);
h_stiff = zeros(num_h_modes,num_h_modes,num_x);
h_force = zeros(num_h_modes,num_x);
for iX = 1:num_x
    r_ddot_i = r_ddot(:,iX);
    r_dot_i = r_dot(:,iX);
    r_transformed_i = r_transformed(:,iX);


    for iTerm = 1:num_coeffs
        r_power_products(iTerm,1) = prod(r_transformed_i.^input_order(:,iTerm));
    end

    r_force = Eom_Input.Reduced_Force_Data.coeffs*r_power_products(1:num_reduced_force_coeffs,:);
    h_stiffness = tensorprod(Eom_Input.H_Stiffness_Data.coeffs,r_power_products(1:num_h_stiffness_coeffs,:),3,1);

    r_products_theta = r_power_products(1:num_coupling_coeffs,:);
    r_dr_products_theta = r_products_theta(Eom_Input.Condensed_Disp_Data.diff_mapping{1,1}).*Eom_Input.Condensed_Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_theta = r_products_theta(Eom_Input.Condensed_Disp_Data.diff_mapping{1,2}).*Eom_Input.Condensed_Disp_Data.diff_scale_factor{1,2};

    r_products_G = r_power_products(1:num_h_coupling_coeffs,:);
    r_dr_products_G = r_products_G(Eom_Input.H_Coupling_Gradient.diff_mapping{1,1}).*Eom_Input.H_Coupling_Gradient.diff_scale_factor{1,1};
    r_dr2_products_G = r_products_G(Eom_Input.H_Coupling_Gradient.diff_mapping{1,2}).*Eom_Input.H_Coupling_Gradient.diff_scale_factor{1,2};

    %%% Inertia
    pre_G_prod = squeeze(tensorprod(r_products_G',G_beta_bar,2,1));
    G_G_prod = tensorprod(pre_G_prod,r_products_G,d2_dims,1);

    h_inertia(:,:,iX) = I_h + G_G_prod;

    %%% Convection
    G_Gdr_prod = tensorprod(pre_G_prod,r_dr_products_G,d2_dims,1);
    h_conv(:,:,iX) = 2*tensorprod(G_Gdr_prod,r_dot_i,3,1);
    
    %%% Stiffness
    stiff_1 = tensorprod(G_Gdr_prod,r_ddot_i,3,1);
    
    G_Gdr2_prod = tensorprod(pre_G_prod,r_dr2_products_G,d2_dims,1);
    pre_stiff_2 = tensorprod(G_Gdr2_prod,r_dot_i,4,1);
    stiff_2 = tensorprod(pre_stiff_2,r_dot_i,3,1);

    h_stiff(:,:,iX) = stiff_1 + stiff_2 + h_stiffness;
    
    %%% Force
    force_r = r_ddot_i + r_force;
    
    L_theta_dr_prod = tensorprod(L_theta_beta,r_dr_products_theta,2,1);
    L_theta_dr2_prod = tensorprod(L_theta_beta,r_dr2_products_theta,2,1);
    force_L_1 = L_theta_dr_prod*r_ddot_i;
    force_L_2 = tensorprod(L_theta_dr2_prod,r_dot_i,3,1)*r_dot_i;
    force_L = force_L_1 + force_L_2;

    pre_G_theta_dr_prod = squeeze(tensorprod(r_products_G',G_theta_beta,2,1));
    G_theta_dr_prod = tensorprod(pre_G_theta_dr_prod,r_dr_products_theta,2,1);
    G_theta_dr2_prod = tensorprod(pre_G_theta_dr_prod,r_dr2_products_theta,2,1);
    pre_force_h_2 = tensorprod(G_theta_dr2_prod,r_dot_i,3,1);
    force_h_1 = G_theta_dr_prod*r_ddot_i;    
    force_h_2 = pre_force_h_2*r_dot_i;
    force_h = force_h_1 + force_h_2;
    
    h_force(:,iX) = -([force_r;force_L] + force_h);
end
end