clear

load("data\exhaust_1567\Static_Data.mat");
%-------------------------------------------------------------------------%
% Static_Data = Static_Data.add_validation_data(129);

Rom = Reduced_System(Static_Data);

ax = Static_Data.plot_force;
Rom.Force_Polynomial.plot_polynomial(ax);

Theta_Poly = Rom.Condensed_Displacement_Polynomial;
% theta_coeffs = Theta_Poly.coefficients;
% r_evecs = Static_Data.Model.reduced_eigenvectors;
% shift_factor = Theta_Poly.shifting_factor;
% scale_factor = Theta_Poly.scaling_factor;
% 
% 
% num_modes = size(r_evecs,2);
% for iMode = 1:num_modes
%     theta_coeffs(1,:) = theta_coeffs(1,:) - r_evecs(:,iMode)'*shift_factor(iMode);
%     theta_coeffs(iMode+1,:) = theta_coeffs(iMode+1,:) - r_evecs(:,iMode)'/scale_factor(iMode);
% end
% Theta_Poly.coefficients = theta_coeffs;

% num_plots = 9;
% plot_ids = randi(8566,[num_plots,1]);
plot_ids = 1804;

ax = Static_Data.plot_condensed_displacement(plot_ids);
Theta_Poly.plot_polynomial(ax,plot_ids)

ax = Static_Data.plot_energy;
Rom.Potential_Polynomial.plot_polynomial(ax)

plot_index = [1323,4];
ax = Static_Data.plot_h_coupling_gradient(plot_index);

% plot_index = randi(1434,[9,1]);
% ax = Static_Data.plot_condensed_displacement(plot_index);
% Rom.Condensed_Displacement_Polynomial.plot_polynomial(ax,plot_index);

% H_Stiffness_Poly = Rom.Low_Frequency_Stiffness_Polynomial;
% H_Coupling_Gradient_Poly = Rom.Low_Frequency_Coupling_Gradient_Polynomial;
% 
% plot_index = [1074,1];
% ax = Static_Data.plot_h_coupling_gradient(plot_index);
% H_Coupling_Gradient_Poly.plot_polynomial(ax,plot_index);
% 
% plot_index = [1,1];
% ax = Static_Data.plot_h_stiffness(plot_index);
% H_Stiffness_Poly.plot_polynomial(ax,plot_index);

