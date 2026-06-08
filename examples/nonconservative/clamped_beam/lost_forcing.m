clear
close all
rom_id = 11001;
Static_Data = load_static_data("clamped_beam_" + string(rom_id));
Rom = Reduced_System(Static_Data);
Model = Rom.Model;
%---
num_dof = Model.num_dof;
y_coords = (1:(num_dof/6)-1)*6 + 2;


displacement = zeros(num_dof,1);
force = zeros(num_dof,1);
force(y_coords(38)) = 0.05;


plot_fe_force(Model,displacement,force);


%--------
vector_rejection = Model.mass*Model.reduced_eigenvectors*Model.reduced_eigenvectors';

remaining_force = vector_rejection*force;
lost_force = force -remaining_force;


% plot_fe_force(Model,displacement,10*remaining_force);
%--------
r = zeros(length(Model.reduced_modes),1);
x = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r);
disp_gradient = Rom.Physical_Displacement_Polynomial.differentiate_polynomial;
x_dr = disp_gradient.evaluate_polynomial(r);
projected_force = x_dr'*force;
lin_projected_force = Model.reduced_eigenvectors'*force;
equivalent_force = Model.mass*Model.reduced_eigenvectors*projected_force;
disp([projected_force,lin_projected_force]);


plot_fe_force(Model,displacement,equivalent_force);
%--------
% num_points = 11;
% r_lim = Rom.reduced_displacement_limits;
% r = linspace(r_lim(1),r_lim(2),num_points);
% disp_gradient = Rom.Physical_Displacement_Polynomial.differentiate_polynomial;


% figure;
% ax = axes;
% hold(ax,"on")
% % ax = plot_fe_force(Model,disp,force,"axes",ax);
% for iPoint = 1:num_points
%     x = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r(:,iPoint));
%     x_dr = disp_gradient.evaluate_polynomial(r(:,iPoint));
%     projected_force = x_dr'*force;
%     lin_projected_force = Model.reduced_eigenvectors'*force;
%     equivalent_force = Model.mass*Model.reduced_eigenvectors*projected_force;
%     disp([projected_force,lin_projected_force]);
% 
% 
%     plot_fe_force(Model,displacement,10*equivalent_force,"axes",ax);
% end
% hold(ax,"off")

