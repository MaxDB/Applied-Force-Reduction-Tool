function x_dot = coco_forced_eom(t,x,force_amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data)
%eom for new force approach
num_x = size(x,2);
num_modes = size(x,1)/2;
num_applied_forces = size(Applied_Force_Data.shape,2);
num_r_modes = num_modes - num_applied_forces;

num_force_coeffs = size(Force_Data.coeffs,2);
num_coupling_coeffs = size(Disp_Data.beta_bar,1);
num_coeffs = size(input_order,1);

disp_span = 1:num_modes;
z = x(disp_span,:);

r_span = 1:num_r_modes;
p_span = (num_r_modes+1):num_modes;


scale_factor = Force_Data.scale_factor;
shift_factor = Force_Data.shift_factor;
%assumes force and coupling from same dataset

z_transformed = scale_factor.*(z + shift_factor);


vel_span = disp_span + num_modes;
z_dot = x(vel_span,:);


x_dot = zeros(2*num_modes,num_x);
x_dot(disp_span,:) = z_dot;

% r_power_products = ones(num_coeffs,1);
frequency = 2*pi./period;
for iX = 1:num_x
    z_i = z_transformed(:,iX);
    t_i = t(iX);
    frequency_i = frequency(iX);
    z_dot_i = z_dot(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_modes
        r_power_products = r_power_products.*z_i(iMode).^input_order(:,iMode);
    end

    r_products_coupling = r_power_products(1:num_coupling_coeffs,:);

    r_dr_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
    r_dr2_products_coupling = r_products_coupling(Disp_Data.diff_mapping{1,2}).*Disp_Data.diff_scale_factor{1,2};
    
    disp_prod = r_dr_products_coupling'*Disp_Data.beta_bar;

    %--
    reduced_restoring_force_z = Force_Data.coeffs*r_power_products(1:num_force_coeffs,:);

    disp_z_prod = r_dr_products_coupling'*Force_Data.disp_z_force_beta;
    restoring_force = disp_z_prod*reduced_restoring_force_z;
    %--
    inertia_term = disp_prod*r_dr_products_coupling;
    %--
    convection_term = disp_prod*(tensorprod(r_dr2_products_coupling,z_dot_i,3,1)*z_dot_i);
    %--
    damping_prod = r_dr_products_coupling'*Damping_Data.damping_beta;
    damping_term = damping_prod*(r_dr_products_coupling)*z_dot_i;
    %--
    disp_amp_prod = r_dr_products_coupling'*Applied_Force_Data.disp_force_beta;
    applied_force = force_amp*disp_amp_prod*sin(frequency_i*t_i);
    %--
    x_dot(vel_span,iX) = -inertia_term\(convection_term+restoring_force + damping_term - applied_force);
end
end