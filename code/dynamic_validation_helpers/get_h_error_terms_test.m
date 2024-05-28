function  [h_inertia,h_conv,h_stiff,h_force] = get_h_error_terms_test_2(r,r_dot,r_ddot,Eom_Input)
num_x = size(r,2);

r_evec = Eom_Input.r_evecs;
L_evec = Eom_Input.h_evecs;

h_evec = [r_evec,L_evec];

num_h_modes = size(h_evec,2);
I_h = eye(num_h_modes);
mass = Eom_Input.mass;



h_inertia = zeros(num_h_modes,num_h_modes,num_x);
h_conv = zeros(num_h_modes,num_h_modes,num_x);
h_stiff = zeros(num_h_modes,num_h_modes,num_x);
h_force = zeros(num_h_modes,num_x);
for iX = 1:num_x
    r_i = r(:,iX);
    r_dot_i = r_dot(:,iX);
    r_ddot_i = r_ddot(:,iX);

    f_r = Eom_Input.Reduced_Force_Poly.evaluate_polynomial(r_i);
    D = Eom_Input.H_Stiffness_Poly.evaluate_polynomial(r_i);

    theta_dr = Eom_Input.Theta_dr_Poly.evaluate_polynomial(r_i)';
    theta_dr2 = Eom_Input.Theta_dr2_Poly.evaluate_polynomial(r_i)';
    

    G = Eom_Input.H_Coupling_Poly.evaluate_polynomial(r_i)';
    G_dr = Eom_Input.H_Coupling_dr_Poly.evaluate_polynomial(r_i)';
    G_dr2 = Eom_Input.H_Coupling_dr2_Poly.evaluate_polynomial(r_i)';


    %%% Inertia
    h_inertia(:,:,iX) = I_h + G'*mass*G;

    %%% Convection
    h_conv(:,:,iX) = 2*G'*mass*(G_dr*r_dot_i);
    
    %%% Stiffness
    h_stiff(:,:,iX) = + G'*mass*(G_dr*r_ddot_i + (G_dr2*r_dot_i)*r_dot_i) + D;
    
    %%% Force
    force_1 = r_ddot_i + f_r;
    force_2 = L_evec'*mass*(theta_dr*r_ddot_i + (theta_dr2*r_dot_i)*r_dot_i);
    force_3 = G'*mass*(theta_dr*r_ddot_i + (theta_dr2*r_dot_i)*r_dot_i);
    
    h_force(:,iX) = -([force_1;force_2] + force_3);
end
end