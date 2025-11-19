function Static_Data = minimum_stiffness_degree(Static_Data)
MAX_DEGREE = 15;

Validation_Opts = Static_Data.Validation_Options;
max_fitting_error = Validation_Opts.maximum_fitting_error;

initial_degree = Validation_Opts.minimum_degree;

V = Static_Data.potential_energy;
validated_points = V < Static_Data.Model.energy_limit;

r = Static_Data.reduced_displacement(:,validated_points);
stiffness_data = Static_Data.tangent_stiffness.nonzero_data(:,validated_points);

stiffness_degree = initial_degree-1;

%find minimum degree

while stiffness_degree  < MAX_DEGREE
    rom = Reduced_System(Static_Data,"degree",[1,1,stiffness_degree]);
    modeled_outputs = rom.Tangent_Stiffness_Polynomial.Nonzero_Polynomials.modeled_outputs;
    modeled_outputs_indicies = find(modeled_outputs);
    
    tangent_stiffness_rom = rom.Tangent_Stiffness_Polynomial.evaluate_polynomial(r,modeled_outputs_indicies);
    stiffness_error = coeff_of_determination(stiffness_data(modeled_outputs,:),tangent_stiffness_rom);
    [max_stiffness_error,max_index] = max(stiffness_error);
    max_stiffness_index = modeled_outputs_indicies(max_index); %#ok<NASGU>
    if max_stiffness_error < max_fitting_error
        break
    end

    stiffness_degree = stiffness_degree + 1;
end
if stiffness_degree >= MAX_DEGREE
    error("Tangent stiffness polynomial requires too high a degree")
end
Static_Data.validated_degree(3) = stiffness_degree;
end