function Validation_Eq_Data = time_parametrise_h_terms(h_terms,t,r,r_dot,r_ddot,period)


[h_mass,h_conv,h_stiff,h_force] = h_terms(t,r,r_dot,r_ddot,period); 

%--
% Validation_Eq_Data.h_mass_spline = spline(t,h_mass);
% Validation_Eq_Data.h_conv_spline = spline(t,h_conv);
% Validation_Eq_Data.h_stiff_spline = spline(t,h_stiff);
% Validation_Eq_Data.h_force_spline = spline(t,h_force);

h_force_mat = permute(h_force,[1,3,2]);

h_stiff_prod = -pagemldivide(h_mass,h_stiff);
h_conv_prod = -pagemldivide(h_mass,h_conv);
h_force_prod = permute(pagemldivide(h_mass,h_force_mat),[1,3,2]);

Validation_Eq_Data.h_stiff_prod_spline = spline(t,h_stiff_prod);
Validation_Eq_Data.h_conv_prod_spline = spline(t,h_conv_prod);
Validation_Eq_Data.h_force_prod_spline = spline(t,h_force_prod);

end