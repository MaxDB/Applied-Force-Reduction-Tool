function max_norm_error = verify_time_parametrised_validation(t,Validation_Eq_Data,Validation_Eq_Data_Verification,Verification_Options)

prob_mat_1 = get_validation_equation(t,Validation_Eq_Data);
prob_mat_2 = get_validation_equation(t,Validation_Eq_Data_Verification);

max_prob_mat = max(max(abs(prob_mat_1),[],3),[],1);
prob_dims = size(prob_mat_1,2);
z_test = eye(prob_dims)./max_prob_mat;

z_dot_1 = pagemtimes(prob_mat_1,z_test);
z_dot_2 = pagemtimes(prob_mat_2,z_test);

ignore_index = z_dot_1 < 1e-3;
z_dot_1(ignore_index) = 0;
z_dot_2(ignore_index) = 0;

all_error = abs((z_dot_1 - z_dot_2)./z_dot_1);
[max_time_error,error_index] = max(all_error,[],3);

max_mode_error = max(max_time_error,[],2);

max_norm_error = max(max_mode_error)/Verification_Options.maximum_interpolation_error(4);

end

function prob_mat = get_validation_equation(t,Validation_Eq_Data)
%--
% h_stiff_prod = ppval(Validation_Eq_Data.h_stiff_prod_spline,t);
% h_conv_prod = ppval(Validation_Eq_Data.h_conv_prod_spline,t);
% 
% prob_dim = size(h_conv_prod,1);
% num_time_points = size(t,2);
% 
% disp_span = 1:prob_dim;
% vel_span = disp_span + prob_dim;
% 
% prob_mat = zeros(prob_dim,2*prob_dim,num_time_points);
% 
% prob_mat(disp_span,disp_span,:) = h_stiff_prod;
% prob_mat(disp_span,vel_span,:) = h_conv_prod;

prob_mat = ppval(Validation_Eq_Data.h_stiff_prod_spline,t);

end

