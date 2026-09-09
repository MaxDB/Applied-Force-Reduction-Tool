function Validated_Trajectory = solve_trajectory_validation(Solution,Validation_Rom,Verification_Rom, Validation_Settings)

num_r_modes = Validation_Rom.get_reduced_dimension;
disp_span = 1:num_r_modes;
vel_span = disp_span + num_r_modes;


Validation_Opts = read_default_options("validation");
Validation_Opts.get_stability = 0;
Validation_Opts.save_orbit = 0;
Validation_Opts.validation_algorithm = "h_ivp";

%%% Set up h-problem
[h_terms,reduced_eom,~,Validation_Analysis_Inputs] = set_up_validation_problem(Validation_Rom,Validation_Opts,Solution,Validation_Settings);
h_terms_verification = set_up_validation_problem(Verification_Rom,Validation_Opts,Solution,Validation_Settings);

%-
forcing_period = 2*pi/Solution.Force_Data.frequency;

t0 = Solution.t;
r = Solution.r;
r_dot = Solution.r_dot;
z_r = reduced_eom(t0,[r;r_dot],forcing_period);
r_ddot = z_r(vel_span,:);

%-
validation_equation_time_start = tic;
%-
Validation_Eq_Data = time_parametrise_h_terms(h_terms,t0,r,r_dot,r_ddot,forcing_period);
Validation_Eq_Data_Verification = time_parametrise_h_terms(h_terms_verification,t0,r,r_dot,r_ddot,forcing_period);
%-
verification_time_start = tic;
error = verify_time_parametrised_validation(t0,Validation_Eq_Data,Validation_Eq_Data_Verification,Validation_Settings.Verification_Options);
verification_time = toc(verification_time_start);
log_message = sprintf("Validation polynomials verification error: %.1f in %.1f seconds" ,error,verification_time);
logger(log_message,3)
%-
validation_equation = @(t,z) get_validation_equation(t,z,Validation_Eq_Data);

%--
x0 = Validation_Rom.expand(r(:,1));
x_dot0 = Validation_Rom.expand_velocity(r(:,1),r_dot(:,1));

r_evecs = Validation_Rom.Model.reduced_eigenvectors;
l_evecs = Validation_Rom.Model.low_frequency_eigenvectors;
if class(r_evecs) == "Large_Matrix_Pointer"
    r_evecs = r_evecs.load();
end
if class(l_evecs) == "Large_Matrix_Pointer"
    l_evecs = l_evecs.load();
end

l_mode_map = ismember(Validation_Rom.Model.low_frequency_modes,Validation_Settings.L_modes);
h_evecs = [r_evecs,l_evecs(:,l_mode_map)];
h_transform = h_evecs'*Validation_Rom.Model.mass;

h0 = h_transform*x0;
h_dot0 = h_transform*x_dot0;
z0 = [h0;h_dot0];
%-
validation_equation_time = toc(validation_equation_time_start);
log_message = sprintf("Validation equations generated in %.1f seconds" ,validation_equation_time);
logger(log_message,2)
%---
validation_sim_time_start = tic;

Validation_Sol = ode45(validation_equation,t0,z0,Solution.ode_options);


num_h_modes = size(h0,1);
disp_span = 1:num_h_modes;
vel_span = disp_span + num_h_modes;

Validated_Trajectory.t = Validation_Sol.x;
Validated_Trajectory.h = Validation_Sol.y(disp_span,:);
Validated_Trajectory.h_dot = Validation_Sol.y(vel_span,:);

r_spline = spline(t0,r);
r_dot_spline = spline(t0,r_dot);


Validated_Trajectory.r = ppval(r_spline,Validation_Sol.x);
Validated_Trajectory.r_dot = ppval(r_dot_spline,Validation_Sol.x);


validation_sim_time = toc(validation_sim_time_start);
log_message = sprintf("Validation trajectory simulated in %.1f seconds" ,validation_sim_time);
logger(log_message,2)
%------
% Analysis
validation_analysis_time_start = tic;
r = Validated_Trajectory.r;
r_dot = Validated_Trajectory.r_dot;
h = Validated_Trajectory.h;
h_dot = Validated_Trajectory.h_dot;

[ke_tilde,ke_hat] = h_kinetic_energy(r,r_dot,h,h_dot,Validation_Analysis_Inputs);

r_force = Validation_Analysis_Inputs.Force_Poly.evaluate_polynomial(r);

num_points = size(r,2);
h_potential = zeros(1,num_points);
for iPoint = 1:num_points
    h_i = h(:,iPoint);
    r_force_hat_i = [r_force(:,iPoint);zeros(num_h_modes-num_r_modes,1)];
    r_force_gradient_i = Validation_Analysis_Inputs.H_Force_Poly.evaluate_polynomial(r(:,iPoint));
    if size(r_force_gradient_i,1) == num_h_modes
        h_potential(iPoint) = (r_force_hat_i + r_force_gradient_i*h_i)'*h_i;
    else
        x_i = Validation_Analysis_Inputs.H_Stiffness_Poly.evaluate_polynomial(r(:,iPoint))*h_i;
        h_potential(iPoint) = r_force_hat_i'*h_i + h_i'*r_force_gradient_i'*x_i;
    end
end

potential_tilde = Validation_Analysis_Inputs.Potential_Poly.evaluate_polynomial(r);
potential_hat = h_potential + potential_tilde;

h_energy = potential_hat + ke_hat;
r_energy = potential_tilde + ke_tilde;

Validated_Trajectory.r_energy = r_energy;
Validated_Trajectory.h_energy = h_energy;

validation_analysis_time = toc(validation_analysis_time_start);
log_message = sprintf("Validation trajectory analysed in %.1f seconds" ,validation_analysis_time);
logger(log_message,2)
end

function z_dot = get_validation_equation(t,z,Validation_Eq_Data)
% h_mass = ppval(Validation_Eq_Data.h_mass_spline,t);
% h_conv = ppval(Validation_Eq_Data.h_conv_spline,t);
% h_stiff = ppval(Validation_Eq_Data.h_stiff_spline,t);
% h_force = ppval(Validation_Eq_Data.h_force_spline,t);
% 
% 
% prob_dim = size(h_force,1);
% disp_span = 1:prob_dim;
% vel_span = disp_span + prob_dim;
% 
% prob_mat = zeros(2*prob_dim);
% prob_mat(disp_span,vel_span) = eye(prob_dim);
% prob_mat(vel_span,disp_span) = -h_mass\h_stiff;
% prob_mat(vel_span,vel_span) = - h_mass\h_conv;
% 
% prob_vec = zeros(2*prob_dim,1);
% prob_vec(vel_span) = h_mass\h_force;
% z_dot = prob_mat*z + prob_vec;

%--
h_stiff_prod = ppval(Validation_Eq_Data.h_stiff_prod_spline,t);
h_conv_prod = ppval(Validation_Eq_Data.h_conv_prod_spline,t);
h_force_prod = ppval(Validation_Eq_Data.h_force_prod_spline,t);

prob_dim = size(h_force_prod,1);
disp_span = 1:prob_dim;
vel_span = disp_span + prob_dim;

prob_mat = zeros(2*prob_dim);
prob_mat(disp_span,vel_span) = eye(prob_dim);
prob_mat(vel_span,disp_span) = h_stiff_prod;
prob_mat(vel_span,vel_span) = h_conv_prod;

prob_vec = zeros(2*prob_dim,1);
prob_vec(vel_span) = h_force_prod;
z_dot = prob_mat*z + prob_vec;

end



