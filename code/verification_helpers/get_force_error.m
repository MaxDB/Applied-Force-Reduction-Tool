function error = get_force_error(disp,Rom_One,Rom_Two)
    potential_one = Rom_One.Potential_Polynomial.evaluate_polynomial(disp);
    potential_two = Rom_Two.Potential_Polynomial.evaluate_polynomial(disp);
    % 
    % error = (2*abs(potential_one - potential_two)./(abs(potential_one) + abs(potential_two)));
    error = abs(potential_one - potential_two)./abs(potential_one);
end