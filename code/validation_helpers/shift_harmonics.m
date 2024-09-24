function shifted_coeffs = shift_harmonics(coeffs,time_shift,frequency)
    a1 = coeffs(2);
    b1 = coeffs(3);

    % phi = frequency*time_shift - pi/2;
    % if a1 == 0
    %     a_shift = b1/sqrt(1+tan(phi)^2);
    %     b_shift = -a_shift*tan(phi);
    %     shifted_coeffs = [coeffs(1),a_shift,b_shift];
    % end

    phi = frequency*time_shift;
    if a1 == 0
        a_shift = b1*sin(phi);
        b_shift = b1*cos(phi);
        shifted_coeffs = [coeffs(1),a_shift,b_shift];
    end

    % t = linspace(0,(2*pi/frequency));
    % x_0 = b1*sin(frequency*t);
    % x_1 = b1*sin(frequency*(t+time_shift));
    % x_2 = a_shift*cos(frequency*t) + b_shift*sin(frequency*t);
end