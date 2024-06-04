clear

load("data\exhaust_157\Static_Data.mat");
%-------------------------------------------------------------------------%
% Static_Data = Static_Data.add_validation_data(129);

Rom = Reduced_System(Static_Data);


% ax = Static_Data.plot_condensed_displacement(plot_ids);
% Theta_Poly.plot_polynomial(ax,plot_ids)
% 
% ax = Static_Data.plot_energy;
% Rom.Potential_Polynomial.plot_polynomial(ax)
% 
% plot_index = [1323,4];
% ax = Static_Data.plot_h_coupling_gradient(plot_index);

% plot_index = randi(1434,[9,1]);
% ax = Static_Data.plot_condensed_displacement(plot_index);
% Rom.Condensed_Displacement_Polynomial.plot_polynomial(ax,plot_index);

H_Stiffness_Poly = Rom.Low_Frequency_Stiffness_Polynomial;
H_Coupling_Gradient_Poly = Rom.Low_Frequency_Coupling_Gradient_Polynomial;
% 
% plot_index = [1074,1];
% ax = Static_Data.plot_h_coupling_gradient(plot_index);
% H_Coupling_Gradient_Poly.plot_polynomial(ax,plot_index);
% 
plot_index = [2,2];
ax = Static_Data.plot_h_stiffness(plot_index);
H_Stiffness_Poly.plot_polynomial(ax,plot_index);

