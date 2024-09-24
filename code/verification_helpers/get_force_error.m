function error = get_force_error(disp,Rom_One,Rom_Two)
    potential_one = Rom_One.Potential_Polynomial.evaluate_polynomial(disp);
    potential_two = Rom_Two.Potential_Polynomial.evaluate_polynomial(disp);
    
    error = abs(2*(potential_one - potential_two)./(potential_one + potential_two));
end