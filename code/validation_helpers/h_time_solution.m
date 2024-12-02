function h_frequency = h_time_solution(h_terms,t0,omega,max_harmonics)

[h_inertia,h_conv,h_stiff,h_force] = h_terms{:};
num_h_modes = size(h_force,1);
I_L = eye(num_h_modes);

num_time_points = length(t0);
num_coefficients = 2*max_harmonics+1;
B = zeros(num_time_points*num_h_modes,num_coefficients*num_h_modes);

harmonic_span = @(index) ((index-1)*num_h_modes+1):(index*num_h_modes);
for iTime = 1:num_time_points
    A_j = zeros(num_h_modes,num_coefficients*num_h_modes);
    A_dot_j = zeros(num_h_modes,num_coefficients*num_h_modes);
    A_ddot_j = zeros(num_h_modes,num_coefficients*num_h_modes);

    %constant terms
    A_j(:,harmonic_span(1)) = I_L;


    for nHarmonic = 1:max_harmonics
        c_nj = I_L*cos(nHarmonic*omega*t0(iTime));
        s_nj = I_L*sin(nHarmonic*omega*t0(iTime));

        term_1 = c_nj;
        term_2 = s_nj;
        A_j(:,harmonic_span(1+nHarmonic)) = term_1;
        A_j(:,harmonic_span(1+max_harmonics+nHarmonic)) = term_2;

        term_1 = -(nHarmonic*omega)*s_nj;
        term_2 = (nHarmonic*omega)*c_nj;
        A_dot_j(:,harmonic_span(1+nHarmonic)) = term_1;
        A_dot_j(:,harmonic_span(1+max_harmonics+nHarmonic)) = term_2;

        term_1 = -(nHarmonic*omega)^2*c_nj;
        term_2 = -(nHarmonic*omega)^2*s_nj;
        A_ddot_j(:,harmonic_span(1+nHarmonic)) = term_1;
        A_ddot_j(:,harmonic_span(1+max_harmonics+nHarmonic)) = term_2;
    end

    B_j = h_inertia(:,:,iTime)*A_ddot_j + h_conv(:,:,iTime)*A_dot_j + h_stiff(:,:,iTime)*A_j;

    B(harmonic_span(iTime),:) = B_j;
end
h_frequency_coefficients = lsqminnorm(B,reshape(h_force,num_time_points*num_h_modes,1));

h_frequency = zeros(num_h_modes,num_coefficients);
for iCoeff = 1:num_coefficients
    h_frequency(:,iCoeff) = h_frequency_coefficients(harmonic_span(iCoeff));
end

end
