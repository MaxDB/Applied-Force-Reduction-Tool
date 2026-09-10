function kinetic_energy = r_kinetic_energy(r,r_dot,Eom_Input)
num_x = size(r,2);
num_r_modes = size(r,1);

num_disp_coeffs = size(Eom_Input.Disp_Data.beta_bar,1);


scale_factor = Eom_Input.Disp_Data.scale_factor;
shift_factor = Eom_Input.Disp_Data.shift_factor;

r_transformed = scale_factor.*(r + shift_factor);

r_disp_beta_bar = Eom_Input.Disp_Data.beta_bar;

input_order = Eom_Input.input_order;
input_order = input_order(1:num_disp_coeffs,:);

kinetic_energy = zeros(1,num_x);
for iX = 1:num_x
    r_dot_i = r_dot(:,iX);
    r_transformed_i = r_transformed(:,iX);

    r_power_products = ones(num_disp_coeffs,1);
    for iMode = 1:num_r_modes
        r_power_products = r_power_products.*r_transformed_i(iMode).^input_order(:,iMode);
    end


    r_products_disp = r_power_products(1:num_disp_coeffs,:);
    r_dr_products_disp = r_products_disp(Eom_Input.Disp_Data.diff_mapping{1,1}).*Eom_Input.Disp_Data.diff_scale_factor{1,1};

 
    %--
    r_dr_r_dot_prod = r_dr_products_disp*r_dot_i;

    %--
    r_disp_r_disp = r_dr_r_dot_prod'*r_disp_beta_bar*r_dr_r_dot_prod;

    kinetic_energy(:,iX) = 0.5*r_disp_r_disp;
end