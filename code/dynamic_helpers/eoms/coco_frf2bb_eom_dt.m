function x_dot_dt = coco_frf2bb_eom_dt(t,x,epsilon,force_amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data)
num_x = size(x,2);
num_modes = size(x,1)/2;

num_coupling_coeffs = size(Disp_Data.beta_bar,1);
num_coeffs = size(input_order,1);

disp_span = 1:num_modes;
r = x(disp_span,:);
scale_factor = Force_Data.scale_factor;
shift_factor = Force_Data.shift_factor;
%assumes force and coupling from same dataset

r_transformed = scale_factor.*(r - shift_factor);


vel_span = disp_span + num_modes;


% r_power_products = ones(num_coeffs,1);
force_type = Applied_Force_Data.type;
switch force_type
    case {"modal","point force"}
        force_shape_dt = Applied_Force_Data.shape_dt(t,force_amp,period);
end

x_dot_dt = zeros(2*num_modes,num_x);
for iX = 1:num_x
    r_i = r_transformed(:,iX);
    epsilon_i = epsilon(iX);

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
    switch force_type
        case "modal"
            applied_force_dt = force_shape_dt(:,iX);
        case "point force"
            amplitude_shape = r_dr_products_coupling'*Applied_Force_Data.disp_force_beta;
            applied_force_dt = amplitude_shape*force_shape_dt(:,iX);
    end
    %--
    nonconservative_term = epsilon_i*applied_force_dt;
    %--
    x_dot_dt(vel_span,iX) = inertia_term\nonconservative_term;
end
end