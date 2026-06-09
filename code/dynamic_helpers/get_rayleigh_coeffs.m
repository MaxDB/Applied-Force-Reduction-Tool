function damping_coeffs = get_rayleigh_coeffs(Model,target_damping,target_modes)
max_mode = max(target_modes);
[~,evals] = eigs(Model.stiffness,Model.mass,max_mode,"smallestabs");

natural_freq = sqrt(diag(evals));
if any(target_modes == 0)
    natural_freq = natural_freq(target_modes ~= 0);
    alpha = 2*natural_freq*target_damping;
    beta = 2*target_damping/natural_freq;

    damping_coeffs = [alpha,beta];
    return
end

natural_freq = natural_freq(target_modes);
A = 2*natural_freq.*target_damping;
B = [ones(size(natural_freq)),natural_freq.^2];
damping_coeffs = B\A;

damping_coeffs = damping_coeffs';

end