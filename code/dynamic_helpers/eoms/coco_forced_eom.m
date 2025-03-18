function x_dot = coco_forced_eom(t,x,force_amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data)
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

x_dot = zeros(2*num_modes,num_x);
x_dot(disp_span,:) = r_dot;

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
    restoring_force = Force_Data.coeffs*r_power_products(1:num_force_coeffs,:);

    r_dr_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,2}).*Disp_Data.diff_scale_factor{1,2};
    
    disp_prod = r_dr_products_coupling'*Disp_Data.beta_bar; 
    %--
    inertia_term = disp_prod*r_dr_products_coupling;
    %--
    convection_term = disp_prod*(tensorprod(r_dr2_products_coupling,r_dot_i,3,1)*r_dot_i);
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
    x_dot(vel_span,iX) = -inertia_term\(convection_term+restoring_force + damping_term - applied_force);
end
end