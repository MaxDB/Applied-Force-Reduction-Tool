function [data,E] = coco_energy(prob,data,u,Potential_Polynomial,Disp_Data,input_index,type)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
r_dot = xcn((num_modes+1):end,:);   %velocity

scale_factor = Disp_Data.scale_factor;
shift_factor = Disp_Data.shift_factor;
%assumes force and coupling from same dataset
r_transformed = scale_factor.*(r + shift_factor);

num_coeffs = size(Disp_Data.beta_bar,2);


switch type
    case "free"
        r_power_products = ones(num_coeffs,1);
        for iMode = 1:num_modes
            r_power_products = r_power_products.*r_transformed(iMode,1).^input_index(1:num_coeffs,iMode);
        end

        r_products = r_power_products(1:num_coeffs,:);
        r_dr_products = r_products(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
        ke_prod = r_dr_products*r_dot(:,1);
        kinetic_energy = 0.5*ke_prod'*Disp_Data.beta_bar*ke_prod;


        potential_energy = Potential_Polynomial.evaluate_polynomial(r(:,1));
        E = kinetic_energy+potential_energy;

    case "forced"
        potential_energy = Potential_Polynomial.evaluate_polynomial(r);
        [~,max_index] = max(potential_energy);
        [~,min_index] = min(potential_energy);

        r_indices = [min_index,max_index];
        num_index = 2;
        kinetic_energy = zeros(1,num_index);
        for iX = 1:num_index
            r_dot_i = r_dot(:,r_indices(iX));

            r_power_products = ones(num_coeffs,1);
            for iMode = 1:num_modes
                r_power_products = r_power_products.*r_transformed(iMode,r_indices(iX)).^input_index(1:num_coeffs,iMode);
            end

            r_products = r_power_products(1:num_coeffs,:);
            r_dr_products = r_products(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
            ke_prod = r_dr_products*r_dot_i;
            kinetic_energy(1,iX) = 0.5*ke_prod'*Disp_Data.beta_bar*ke_prod;
        end

        E = max(potential_energy(1,r_indices) + kinetic_energy);
end



end