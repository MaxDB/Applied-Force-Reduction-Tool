function [data,E] = coco_energy(prob,data,u,Potential_Polynomial,Disp_Data,input_index)
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
r_transformed = scale_factor.*(r - shift_factor);

num_coeffs = size(Disp_Data.beta_bar,2);

% input_order = input_index';
% r_power_products = ones(num_coeffs,1);
% 
% 
% for iTerm = 1:num_coeffs
%     r_power_products(iTerm,1) = prod(r_transformed(:,1).^input_order(:,iTerm));
% end

r_power_products = ones(num_coeffs,1);
for iMode = 1:num_modes
    r_power_products = r_power_products.*r_transformed(iMode,1).^input_index(:,iMode);
end

r_products = r_power_products(1:num_coeffs,:);
r_dr_products = r_products(Disp_Data.diff_mapping{1,1}).*Disp_Data.diff_scale_factor{1,1};
kinetic_energy = 0.5*r_dot(:,1)'*r_dr_products'*Disp_Data.beta_bar*r_dr_products*r_dot(:,1) ...
    +0.5*r_dot(:,1)'*r_dot(:,1);


potential_energy = Potential_Polynomial.evaluate_polynomial(r(:,1));

E = kinetic_energy+potential_energy;

end