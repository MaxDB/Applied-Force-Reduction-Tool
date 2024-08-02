function perturbation_scale_factor = select_perturbation_scale_factor(Model) 
r_evals = Model.reduced_eigenvalues;
L_evals = Model.low_frequency_eigenvalues;
h_evals = [r_evals;L_evals];

v_lim = Model.energy_limit;

perturbation_scale_factor = sqrt(2*h_evals*v_lim)';

end