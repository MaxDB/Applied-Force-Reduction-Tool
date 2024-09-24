function error = coeff_of_determination(y,f,num_variables)
sample_size = size(y,2);

ss_res = sum((y-f).^2,2);
y_mean = mean(y,2);
ss_tot = sum((y-y_mean).^2,2);

error = ss_res./ss_tot;
% adjusted_error = error*(sample_size-1)/(sample_size - num_variables - 1);
end