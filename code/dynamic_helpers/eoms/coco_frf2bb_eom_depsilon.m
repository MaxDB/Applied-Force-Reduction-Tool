function x_dot_depsilon = coco_frf2bb_eom_depsilon(t,x,epsilon,force_amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;
num_force_coeffs = size(Force_Data.coeffs,2);
num_disp_coeffs = size(Disp_Data.beta_bar,1);
num_coeffs = size(input_order,1);

disp_span = 1:num_modes;
r = x(disp_span,:);



scale_factor = Force_Data.scale_factor;
shift_factor = Force_Data.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r + shift_factor);


vel_span = disp_span + num_modes;
r_dot = x(vel_span,:);


x_dot_depsilon = zeros(2*num_modes,1,num_x);

% r_power_products = ones(num_coeffs,1);
frequency = 2*pi./period;
for iX = 1:num_x
    r_i = r_transformed(:,iX);
    t_i = t(iX);
    r_dot_i = r_dot(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
    end

    r_products_disp = r_power_products(1:num_disp_coeffs,:);

    r_dr_products_disp = r_products_disp(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};

    disp_prod = r_dr_products_disp'*Disp_Data.beta_bar;
    %--
  
    %--
    inertia_term = disp_prod*r_dr_products_disp;
    %--
    %--
    switch Damping_Data.type
        case "matrix"
            damping_prod = r_dr_products_disp'*Damping_Data.damping_beta;
            damping_term = damping_prod*(r_dr_products_disp)*r_dot_i;
        case "nonlinear_rayleigh"
             r_products_force = r_power_products(1:num_force_coeffs,:);
            r_dr_products_force = r_products_force(Force_Data.diff_mapping{1,1}).*Force_Data.diff_scale_factor{1,1};
            stiffness = Force_Data.coeffs*r_dr_products_force;
            stiffness_damping_prod = r_dr_products_disp'*Damping_Data.disp_r_mode_beta;
            stiffness_damping_term = stiffness_damping_prod*stiffness;

            mass_damping_term = inertia_term;

            damping_term = (Damping_Data.coeffs(1)*mass_damping_term + Damping_Data.coeffs(2)*stiffness_damping_term)*r_dot_i;
    end
    %--
    disp_amp_prod = r_dr_products_disp'*Applied_Force_Data.disp_force_beta;
    applied_force = force_amp*disp_amp_prod*sin(frequency*t_i);
    %--
    x_dot_depsilon(vel_span,1,iX) = -inertia_term\(damping_term - applied_force);
end

end