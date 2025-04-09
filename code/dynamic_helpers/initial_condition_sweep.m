function initial_condition_sweep(Rom)
SWEEP_RESOLUTION = 10;
omega_n = sqrt(Rom.Model.reduced_eigenvalues(1));
period = 2*pi/omega_n;
max_period = 2*period;
t_span = [0,max_period];

eom = Rom.get_equation_of_motion();
input_limits = Rom.Force_Polynomial.input_limit;
r_lim = [min(input_limits,[],2),max(input_limits,[],2)];

num_inputs = size(r_lim,1);

possible_ics = zeros(num_inputs,SWEEP_RESOLUTION);
for iInput = 1:num_inputs
    possible_ics(iInput,:) = linspace(r_lim(iInput,1),r_lim(iInput,2),SWEEP_RESOLUTION);
end

index = ones(num_inputs,1);
num_ics = SWEEP_RESOLUTION^num_inputs;

for iIC = 1:num_ics
    ic = zeros(2*num_inputs,1);
    for iInput = 1:num_inputs
        ic(iInput) = possible_ics(iInput,index(iInput));
    end

    [t,z] = ode45(@(t,z) eom(t,z),t_span,ic);

    index = increment_index(index,SWEEP_RESOLUTION);
end


end

function index = increment_index(index,max_index)
index_size = size(index,1);
counter = index_size;
while index(counter) == max_index
    index(counter) = 1;
    counter = counter - 1;
    if counter == 0
        return
    end
end
index(counter) = index(counter) + 1; 
end