function x_dot_depsilon = coco_frf2bb_eom_depsilon(t,x,epsilon,force_amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

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

x_dot_depsilon = zeros(2*num_modes,1,num_x);

% r_power_products = ones(num_coeffs,1);
force_type = Applied_Force_Data.type;
switch force_type
    case {"modal","point force"}
        force_shape = Applied_Force_Data.shape(t,force_amp,period);
end

for iX = 1:num_x
    r_i = r_transformed(:,iX);
    r_dot_i = r_dot(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
    end

    r_products_coupling = r_power_products(1:num_coupling_coeffs,:);

    r_dr_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
    
    disp_prod = r_dr_products_coupling'*Disp_Data.beta_bar; 
    %--
    inertia_term = disp_prod*r_dr_products_coupling;
    %--
    %--
    damping_prod = r_dr_products_coupling'*Damping_Data.damping_beta;
    damping_term = damping_prod*(r_dr_products_coupling)*r_dot_i;
    %--
    switch force_type
        case "modal"
            applied_force = force_shape(:,iX);
        case "point force"
            amplitude_shape = r_dr_products_coupling'*Applied_Force_Data.disp_force_beta;
            applied_force = amplitude_shape*force_shape(:,iX);
    end
    %--
    nonconservative_term_depsilon = damping_term - applied_force;
    %--
    x_dot_depsilon(vel_span,1,iX) = -inertia_term\nonconservative_term_depsilon;
end
end