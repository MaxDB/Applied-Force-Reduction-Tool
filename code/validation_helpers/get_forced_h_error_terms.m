function  [h_inertia,h_conv,h_stiff,h_force] = get_forced_h_error_terms(t,r,r_dot,r_ddot,amp,period,Eom_Input,epsilon)
num_x = size(r,2);
num_r_modes = size(r,1);

num_reduced_force_coeffs = size(Eom_Input.Reduced_Force_Data.coeffs,2);
num_disp_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp_r_disp,3);
num_h_force_grad_coeffs = size(Eom_Input.H_Force_Data.coeffs,3);
num_h_disp_grad_coeffs = size(Eom_Input.Beta_Bar_Data.h_disp,1);
num_coeffs = size(Eom_Input.input_order,1);

num_h_modes = size(Eom_Input.H_Force_Data.coeffs,2);
num_L_modes = num_h_modes - num_r_modes;

scale_factor = Eom_Input.scale_factor;
shift_factor = Eom_Input.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r + shift_factor);

beta_h_disp = Eom_Input.Beta_Bar_Data.h_disp;
beta_h_disp_r_disp = Eom_Input.Beta_Bar_Data.h_disp_r_disp;

beta_h_disp_h_mode = Eom_Input.Frf_Beta_Bar_Data.h_disp_force_grad;
beta_h_disp_r_mode = Eom_Input.Beta_Bar_Data.h_disp_r_force;

switch Eom_Input.Damping_Data.type
    case "matrix"
        beta_h_disp_damping = Eom_Input.Frf_Beta_Bar_Data.h_disp_damp;
        beta_h_disp_damping_r_disp = Eom_Input.Frf_Beta_Bar_Data.h_disp_damp_r_disp;
    case "nonlinear_rayleigh"
        beta_r_disp_mass_evec_h = Eom_Input.Damping_Beta_Bar.r_disp_mass_evec_h;
        beta_h_disp_mass_evec_h = Eom_Input.Frf_Beta_Bar_Data.h_disp_force_grad;
        beta_r_disp_mass_evec_r = beta_r_disp_mass_evec_h(:,1:num_r_modes);
        beta_h_disp_mass_evec_r = beta_h_disp_mass_evec_h(:,:,1:num_r_modes);

        damping_coeffs = Eom_Input.Damping_Data.coeffs;
end

beta_h_disp_applied_force = Eom_Input.Frf_Beta_Bar_Data.h_disp_applied_force;


%- Apply epsilon
amp = epsilon*amp;
switch Eom_Input.Damping_Data.type
    case "matrix"
        beta_h_disp_damping = epsilon*beta_h_disp_damping;
        beta_h_disp_damping_r_disp = epsilon*beta_h_disp_damping_r_disp;
    case "nonlinear_rayleigh"
        damping_coeffs = epsilon*damping_coeffs;
end
%--


Applied_Force_Data = Eom_Input.Applied_Force_Data;
force_type = Applied_Force_Data.type;
switch force_type
    case {"modal","point force","shape"}
        force_shape = Applied_Force_Data.shape(t,amp,period);
end

input_order = Eom_Input.input_order;

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


    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_r_modes
        r_power_products = r_power_products.*r_transformed_i(iMode).^input_order(:,iMode);
    end

    %---
   
    
    %---
    r_products_force = r_power_products(1:num_reduced_force_coeffs,:);
    r_force = Eom_Input.Reduced_Force_Data.coeffs*r_products_force;

    h_force_grad = tensorprod(Eom_Input.H_Force_Data.coeffs,r_power_products(1:num_h_force_grad_coeffs,:),3,1);
    %--

    r_products_disp = r_power_products(1:num_disp_coeffs,:);
    r_dr_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,1}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_disp = r_products_disp(Eom_Input.Physical_Disp_Data.diff_mapping{1,2}).*Eom_Input.Physical_Disp_Data.diff_scale_factor{1,2};
    
    r_products_grad = r_power_products(1:num_h_disp_grad_coeffs);
    r_dr_products_grad = r_products_grad(Eom_Input.Disp_Grad_Data.diff_mapping{1,1}).*Eom_Input.Disp_Grad_Data.diff_scale_factor{1,1};
    r_dr2_products_grad = r_products_grad(Eom_Input.Disp_Grad_Data.diff_mapping{1,2}).*Eom_Input.Disp_Grad_Data.diff_scale_factor{1,2};
    
    %%% Damping terms

    
    %%% Inertia
    H_beta_prod = squeeze(tensorprod(r_products_grad',beta_h_disp,2,1));
    h_inertia_prod = tensorprod(H_beta_prod,r_products_grad,3,1);
    h_inertia(:,:,iX) = h_inertia_prod;


    %%% Convection
    r_dr_r_dot_prod = r_dr_products_grad*r_dot_i;
    H_r_dr_r_dot_prod = tensorprod(H_beta_prod,r_dr_r_dot_prod,3,1);
    
    h_conv(:,:,iX) = 2*H_r_dr_r_dot_prod;
    
    %%% Stiffness
    r_d2r_r_dot_r_dot_prod = tensorprod(r_dr2_products_grad,r_dot_i,3,1)*r_dot_i;
    r_dr_r_ddot_prod = r_dr_products_grad*r_ddot_i;
    stiff_sum = r_d2r_r_dot_r_dot_prod + r_dr_r_ddot_prod;
    stiff_prod = tensorprod(H_beta_prod,stiff_sum,3,1);


    h_force_grad_pre = tensorprod(beta_h_disp_h_mode,h_force_grad,3,1);
    h_force_grad_projected = squeeze(tensorprod(r_products_grad',h_force_grad_pre,2,1));

    h_stiff(:,:,iX) = stiff_prod + h_force_grad_projected;
    
    %%% Force
    H_beta_disp_prod = squeeze(tensorprod(r_products_grad',beta_h_disp_r_disp,2,1));
    disp_r_d2r_r_dot_r_dot_prod = tensorprod(r_dr2_products_disp,r_dot_i,3,1)*r_dot_i;
    disp_r_dr_r_ddot_prod = r_dr_products_disp*r_ddot_i;
    disp_stiff_sum = disp_r_d2r_r_dot_r_dot_prod + disp_r_dr_r_ddot_prod;
    force_1 = H_beta_disp_prod*disp_stiff_sum;

    %!%
    force_2_pre = tensorprod(beta_h_disp_r_mode,r_force,3,1);
    force_2 = squeeze(tensorprod(r_products_grad',force_2_pre,2,1))';


    conservative_h_force = force_1 + force_2;

    switch force_type
        case "modal"
            applied_reduced_force = force_shape(:,iX);
            applied_force = [applied_reduced_force;zeros(num_L_modes,1)];
        case "point force"
            amplitude_shape = r_products_disp'*Applied_Force_Data.h_disp_force_beta;
            applied_force = amplitude_shape'*force_shape(:,iX);
        case "shape"
            amplitude_shape = r_products_grad'*beta_h_disp_applied_force;
            amplitude_shape = amplitude_shape';
            applied_force = amplitude_shape*force_shape(:,iX);
    end

    h_force(:,iX) = applied_force - conservative_h_force;
    
    %Damping
    switch Eom_Input.Damping_Data.type
        case "matrix"
            pre_H_damping_H_prod = squeeze(tensorprod(r_products_grad',beta_h_disp_damping,2,1));
            H_damping_H = tensorprod(pre_H_damping_H_prod,r_products_grad,d2_dims,1);
            damping_conv = H_damping_H;

            H_damping_Hdr_prod = tensorprod(pre_H_damping_H_prod,r_dr_products_grad,d2_dims,1);
            H_damping_Hdr_r_dot = tensorprod(H_damping_Hdr_prod,r_dot_i,3,1);
            damping_stiff = H_damping_Hdr_r_dot;

            pre_H_damping_disp_prod = squeeze(tensorprod(r_products_grad',beta_h_disp_damping_r_disp,2,1));
            H_damping_disp_prod = pre_H_damping_disp_prod*r_dr_products_disp;
            H_damping_disp_r_dot_prod = H_damping_disp_prod*r_dot_i;
            damping_force = H_damping_disp_r_dot_prod;
        case "nonlinear_rayleigh"
            % beta_r_disp_mass_evec_h = Eom_Input.Damping_Beta_Bar.r_disp_mass_evec_h;
            beta_h_disp_mass_evec_h = Eom_Input.Frf_Beta_Bar_Data.h_disp_force_grad;
            % beta_r_disp_mass_evec_r = beta_r_disp_mass_evec_h(:,1:num_r_modes);
            % beta_h_disp_mass_evec_r = beta_h_disp_mass_evec_h(:,:,1:num_r_modes);
            r_dr_products_force = r_products_force(Eom_Input.Reduced_Force_Data.diff_mapping{1,1}).*Eom_Input.Reduced_Force_Data.diff_scale_factor{1,1};
            reduced_stiffness = zeros(num_h_modes);
            reduced_stiffness(1:num_r_modes,1:num_r_modes) = Eom_Input.Reduced_Force_Data.coeffs*r_dr_products_force;

            h_disp_mass_evec_h = squeeze(tensorprod(r_products_grad',beta_h_disp_mass_evec_h,2,1));
            h_disp_dr_mass_evec_h = tensorprod(r_dr_products_grad',beta_h_disp_mass_evec_h,2,1);
            r_dot_h_disp_dr_mass_evec_h = squeeze(tensorprod(r_dot_i',h_disp_dr_mass_evec_h,2,1));
            r_disp_dr_mass_evec_h = r_dr_products_disp'*beta_r_disp_mass_evec_h;
            r_dot_r_disp_dr_mass_evec_h = r_dot_i'*r_disp_dr_mass_evec_h;

            mass_damping_conv = h_inertia_prod;
            mass_damping_stiff = H_r_dr_r_dot_prod;
            mass_damping_force = H_beta_disp_prod*r_dr_products_disp*r_dot_i;
            
            stiff_damping_prod = h_disp_mass_evec_h*(h_force_grad_projected);
            % stiff_damping_prod = h_disp_mass_evec_h*(reduced_stiffness + h_force_grad_projected);
            stiff_damping_conv = stiff_damping_prod*h_disp_mass_evec_h';
            stiff_damping_stiff = stiff_damping_prod*r_dot_h_disp_dr_mass_evec_h';
            stiff_damping_force = stiff_damping_prod*r_dot_r_disp_dr_mass_evec_h';

            damping_conv = damping_coeffs(1)*mass_damping_conv + damping_coeffs(2)*stiff_damping_conv;
            damping_stiff = damping_coeffs(1)*mass_damping_stiff + damping_coeffs(2)*stiff_damping_stiff;
            damping_force = damping_coeffs(1)*mass_damping_force + damping_coeffs(2)*stiff_damping_force;
    end

    h_force(:,iX) = h_force(:,iX) - damping_force;
    h_stiff(:,:,iX) = h_stiff(:,:,iX)  + damping_stiff;
    h_conv(:,:,iX) = h_conv(:,:,iX) + damping_conv;

end
end