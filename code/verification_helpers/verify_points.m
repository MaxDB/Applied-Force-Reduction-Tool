function [force_error,disp_error] = verify_points(r_points,Rom_One,Rom_Two)
Model = Rom_One.Model;

Disp_Error_Inputs.beta_bar_one = Rom_One.get_beta_bar(Rom_One.Physical_Displacement_Polynomial);
Disp_Error_Inputs.beta_bar_two = Rom_Two.get_beta_bar(Rom_Two.Physical_Displacement_Polynomial);

disp_coeffs_one = Rom_One.Physical_Displacement_Polynomial.coefficients;
disp_coeffs_two = Rom_Two.Physical_Displacement_Polynomial.coefficients;

force_coeffs_one = Rom_One.Force_Polynomial.coefficients;
force_coeffs_two = Rom_Two.Force_Polynomial.coefficients;

force_transform = Model.mass*Model.reduced_eigenvectors;
Disp_Error_Inputs.disp_mode_beta_one = disp_coeffs_one'*force_transform*force_coeffs_one;
Disp_Error_Inputs.disp_mode_beta_two = disp_coeffs_two'*force_transform*force_coeffs_two;

Disp_Error_Inputs.input_order = Rom_Two.get_max_input_order;


Disp_Error_Inputs.Disp_Diff_Data_One = Rom_One.Physical_Displacement_Polynomial.get_diff_data(1);
Disp_Error_Inputs.Disp_Diff_Data_Two = Rom_Two.Physical_Displacement_Polynomial.get_diff_data(1);


restoring_force = Rom_One.Force_Polynomial.evaluate_polynomial(r_points);
force_norm = vecnorm(restoring_force);
force_ratio = restoring_force./force_norm;
force_ratio(abs(force_ratio) < 0.1) = 0;

force_error = get_force_error(r_points,Rom_One,Rom_Two);
disp_error = get_disp_error(r_points,Rom_One,Rom_Two,force_ratio,Disp_Error_Inputs);

end