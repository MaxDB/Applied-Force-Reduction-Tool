function stiffness_error = get_stiffness_error(disp,validation_points,Rom_One,Rom_Two)
num_points = size(disp,2);


num_test_cases = size(validation_points,2);
stiffness_error_all = zeros(num_points,num_test_cases);

for iDisp = 1:num_points
    stiffness_one = Rom_One.Low_Frequency_Stiffness_Polynomial.evaluate_polynomial(disp(:,iDisp));
    stiffness_two = Rom_Two.Low_Frequency_Stiffness_Polynomial.evaluate_polynomial(disp(:,iDisp));
    for iTest = 1:num_test_cases
        h_test = validation_points(:,iTest)*2;
        energy_one = h_test'*stiffness_one*h_test;
        energy_two = h_test'*stiffness_two*h_test;

        stiffness_error_all(iDisp,iTest) = abs(energy_one - energy_two)/abs(energy_one);
    end
end

stiffness_error = max(stiffness_error_all,[],2);
end