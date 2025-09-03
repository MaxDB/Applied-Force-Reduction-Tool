function error = get_disp_error(disp,Rom_One,Rom_Two,force_ratio,Disp_Error_Inputs)
    beta_bar_one = Disp_Error_Inputs.beta_bar_one;
    beta_bar_two = Disp_Error_Inputs.beta_bar_two;
    input_order = Disp_Error_Inputs.input_order;
    Disp_Diff_One = Disp_Error_Inputs.Disp_Diff_Data_One;
    Disp_Diff_Two = Disp_Error_Inputs.Disp_Diff_Data_Two;

    scale_factor = Rom_One.Force_Polynomial.scaling_factor;
    shift_factor = Rom_One.Force_Polynomial.shifting_factor;
    r_transformed = scale_factor.*(disp + shift_factor);

    num_disp_coeffs = [size(beta_bar_one,1),size(beta_bar_two,1)];
    num_force_coeffs = [Rom_One.Force_Polynomial.num_element_coefficients,Rom_Two.Force_Polynomial.num_element_coefficients];
    num_coeffs = max(max(num_disp_coeffs),max(num_force_coeffs));

    num_x = size(r_transformed,2);
    num_modes = size(r_transformed,1);

    error_map = ~ismember(force_ratio,0);
    r_ddot_one = zeros(num_modes,num_x);
    point_error = zeros(num_modes,num_x);
    for iX = 1:num_x
        r_i = r_transformed(:,iX);

        r_power_products = ones(num_coeffs,1);
        for iMode = 1:num_modes
            r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
        end
        
        r_products_disp_one = r_power_products(1:num_disp_coeffs(1),:);
        r_products_disp_two = r_power_products(1:num_disp_coeffs(2),:);

        r_dr_products_disp_one = r_products_disp_one(Disp_Diff_One.diff_mapping{1,1}).*Disp_Diff_One.diff_scale_factor{1,1};
        r_dr_products_disp_two = r_products_disp_two(Disp_Diff_Two.diff_mapping{1,1}).*Disp_Diff_Two.diff_scale_factor{1,1};

        r_products_force_one = r_power_products(1:num_force_coeffs(1),:);
        r_products_force_two = r_power_products(1:num_force_coeffs(2),:);
    
        %----------
        force_one = Rom_One.Force_Polynomial.coefficients*r_products_force_one;
        force_two = Rom_Two.Force_Polynomial.coefficients*r_products_force_two;

        %---------
        inertia_one = r_dr_products_disp_one'*beta_bar_one*r_dr_products_disp_one;
        inertia_two = r_dr_products_disp_two'*beta_bar_two*r_dr_products_disp_two;

        %---------
        r_ddot_one(:,iX) = - inertia_one\force_one;
        r_ddot_two = - inertia_two\force_two;
        % 
        % point_error(:,iX) = 2*abs((r_ddot_one(:,iX) - r_ddot_two)./(abs(r_ddot_one(:,iX)) + abs(r_ddot_two)));

        point_error(:,iX) = abs(r_ddot_one(:,iX) - r_ddot_two)./abs(r_ddot_one(:,iX));
        point_error(~error_map,iX) = 0;
  
    end

    r_ddot_max = max(abs(r_ddot_one),[],2);
    largest_r_ddot = max(r_ddot_max);
    small_acceleration = r_ddot_max < largest_r_ddot/10;
    point_error(small_acceleration,:) = 0;
    error = max(point_error,[],1);
end