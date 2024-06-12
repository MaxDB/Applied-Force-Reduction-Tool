function [h_stiffness,h_stiffness_0,h_coupling_gradient,h_coupling_gradient_0] = parse_perturbation_data(Static_Data,L_modes)
Model = Static_Data.Model;
r_modes = Model.reduced_modes;
num_r_modes = size(r_modes,2);
r_evec = Model.reduced_eigenvectors;

all_L_modes = Model.low_frequency_modes;
[~,L_map] = ismember(L_modes,all_L_modes);
L_evec = Model.low_frequency_eigenvectors(:,L_map);

h_modes = [r_modes,L_modes];
num_h_modes = size(h_modes,2);
h_evec = [r_evec,L_evec];
h_map = [1:num_r_modes,L_map+num_r_modes];


num_dofs = Model.num_dof;
mass = Model.mass;
h_disp_transform = h_evec'*mass;

lambda = Static_Data.perturbation_scale_factor;
F_h = lambda*eye(num_h_modes);

x_perturbation = Static_Data.perturbation_displacement(:,h_map,:);
num_loadcases = size(x_perturbation,3);

%-- setup r perturbations
% rom_degree = Static_Data.validated_degree;
% rom_degree(3:4) = 1;
% Rom = Reduced_System(Static_Data,rom_degree);
% Stiffness_Poly = Rom.Reduced_Stiffness_Polynomial;
% Theta_Tilde_Poly = Rom.Condensed_Displacement_Polynomial;
% Theta_Gradient_Poly = Theta_Tilde_Poly.differentiate_polynomial;
% F_r = lambda*eye(num_r_modes);
% 
% r = Static_Data.reduced_displacement;
%--


h_coupling_gradient = zeros(num_dofs,num_h_modes,num_loadcases);
h_stiffness = zeros(num_h_modes,num_h_modes,num_loadcases);
parfor iLoad = 1:num_loadcases
    % r_i = r(:,iLoad);
    % stiff_i = Stiffness_Poly.evaluate_polynomial(r_i);
    % theta_gradient_i = Theta_Gradient_Poly.evaluate_polynomial(r_i);
    % r_modal_perturbation_i = stiff_i\F_r;
    % r_perturbation_i = (r_evec+theta_gradient_i)*r_modal_perturbation_i;

    % L_perturbation_i = x_perturbation(:,:,iLoad);
    % h_perturbation_i = [r_perturbation_i,L_perturbation_i];
    

    h_perturbation_i = x_perturbation(:,:,iLoad);

    h_disp = h_disp_transform*h_perturbation_i;

    theta_hat = h_perturbation_i - h_evec*h_disp;
    h_coupling_gradient(:,:,iLoad) = theta_hat/h_disp;

    h_stiffness(:,:,iLoad) =  F_h/h_disp;

end

stiffness = Model.stiffness;
h_perturbation_0 = lambda*(stiffness\(h_disp_transform'));
h_0 = h_disp_transform*h_perturbation_0;
h_0 = diag(diag(h_0));
theta_hat_0 = h_perturbation_0 - h_evec*h_0;

h_coupling_gradient_0 = theta_hat_0/h_0;
h_stiffness_0 =  F_h/h_0;
end