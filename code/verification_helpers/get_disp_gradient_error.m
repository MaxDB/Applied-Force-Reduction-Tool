function disp_grad_error = get_disp_gradient_error(disp,validation_points,Rom_One,Rom_Two,force_ratio,Disp_Error_Inputs)

beta_bar_h_disp_one = Disp_Error_Inputs.Beta_Bar_Data_One.h_disp;
beta_bar_h_disp_two = Disp_Error_Inputs.Beta_Bar_Data_Two.h_disp;

input_order = Disp_Error_Inputs.input_order;


scale_factor = Rom_One.Force_Polynomial.scaling_factor;
shift_factor = Rom_One.Force_Polynomial.shifting_factor;
r_transformed = scale_factor.*(disp + shift_factor);

num_force_coeffs = [Rom_One.Force_Polynomial.num_element_coefficients,Rom_Two.Force_Polynomial.num_element_coefficients];
num_h_force_grad_coeffs = [Rom_One.Low_Frequency_Stiffness_Polynomial.num_element_coefficients,Rom_Two.Low_Frequency_Stiffness_Polynomial.num_element_coefficients];
num_h_coupling_grad_coeffs = [size(beta_bar_h_disp_one,1),size(beta_bar_h_disp_two,1)];
num_coeffs = size(input_order,1);


num_x = size(r_transformed,2);
num_r_modes = size(r_transformed,1);
num_h_modes = size(Rom_One.Low_Frequency_Stiffness_Polynomial,1);

num_test_points = size(validation_points,2);
error_map = ~ismember(validation_points,0);

num_points = size(disp,2);
disp_grad_error_all = zeros(num_points,num_test_points);

for iX = 1:num_x
    r_transformed_i = r_transformed(:,iX);

    r_power_products = ones(num_coeffs,1);
    for iMode = 1:num_r_modes
        r_power_products = r_power_products.*r_transformed_i(iMode).^input_order(:,iMode);
    end

    %---
    
    %%% Force
    r_force_one = Rom_One.Force_Polynomial.coefficients*r_power_products(1:num_force_coeffs(1),:);
    r_force_two = Rom_Two.Force_Polynomial.coefficients*r_power_products(1:num_force_coeffs(2),:);

    force_one = zeros(num_h_modes,1);
    force_two = zeros(num_h_modes,1);

    force_one(1:num_r_modes,1) = r_force_one;
    force_two(1:num_r_modes,1) = r_force_two;

    %%% Stiffness
    h_force_grad_one = tensorprod(Rom_One.Low_Frequency_Stiffness_Polynomial.coefficients,r_power_products(1:num_h_force_grad_coeffs(1),:),3,1);
    h_force_grad_two = tensorprod(Rom_Two.Low_Frequency_Stiffness_Polynomial.coefficients,r_power_products(1:num_h_force_grad_coeffs(2),:),3,1);
    
    r_products_disp_one = r_power_products(1:num_h_coupling_grad_coeffs(1),:);
    r_products_disp_two = r_power_products(1:num_h_coupling_grad_coeffs(2),:);
    %%% Inertia
    H_beta_prod_one = squeeze(tensorprod(r_products_disp_one',beta_bar_h_disp_one,2,1));
    H_beta_prod_two = squeeze(tensorprod(r_products_disp_two',beta_bar_h_disp_two,2,1));

    h_inertia_one = tensorprod(H_beta_prod_one,r_products_disp_one,3,1);
    h_inertia_two = tensorprod(H_beta_prod_two,r_products_disp_two,3,1);

    %%%
    pre_acceleration_one = - h_inertia_one\h_force_grad_one;
    pre_acceleration_two = - h_inertia_two\h_force_grad_two;
    
    force_acceleration_one = - h_inertia_one\force_one;
    force_acceleration_two = - h_inertia_two\force_two;
    for iTest = 1:num_test_points
        h_test = validation_points(:,iTest);
        acceleration_one = pre_acceleration_one*h_test + force_acceleration_one;
        acceleration_two = pre_acceleration_two*h_test + force_acceleration_two;
        
        acceleration_error = abs(acceleration_one - acceleration_two)./abs(acceleration_one);
        acceleration_error(~error_map(:,iTest)) = 0;
        
        h_ddot_max = max(abs(acceleration_one),[],2);
        largest_h_ddot = max(h_ddot_max);
        small_acceleration = h_ddot_max < largest_h_ddot/10;
        acceleration_error(small_acceleration) = 0;

        disp_grad_error_all(iX,iTest) = max(acceleration_error);
    end

end

disp_grad_error = max(disp_grad_error_all,[],2);
end