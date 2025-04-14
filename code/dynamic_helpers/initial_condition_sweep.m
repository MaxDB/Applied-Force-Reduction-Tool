function initial_condition_sweep(Rom)
SWEEP_RESOLUTION = 21;
MAX_ERROR = 1e-3;

omega_n = sqrt(Rom.Model.reduced_eigenvalues(1));
period = 2*pi/omega_n;
max_period = 3*period;

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

figure
hold on
low_error_ic = [];
for iIC = 1:num_ics
    ic = zeros(2*num_inputs,1);
    for iInput = 1:num_inputs
        ic(iInput) = possible_ics(iInput,index(iInput));
    end

    [t,z] = ode45(@(t,z) eom(t,z),linspace(0,max_period,50),ic);
    z = z';
    t = t';
    fundamental_frequency = check_fundamental_frequency(t,z');
    periodicity_error = get_periodicity_error(t,z(1:num_inputs,:),fundamental_frequency(1:num_inputs));
    if all(periodicity_error < MAX_ERROR)
        low_error_ic(:,end+1) = ic(1:num_inputs)'; %#ok<AGROW>
    end
    plot3(ic(1),ic(2),min(periodicity_error),"x")
    % plot3(ic(1),ic(2),max(periodicity_error),"x")
    index = increment_index(index,SWEEP_RESOLUTION);
end
hold off

end

function fundamental_frequency = check_fundamental_frequency(t,z)
    z_frequency = fft(z);

    signal_length = size(t,2);
    samping_frequency = 2*pi/diff(t(1:2));
    frequency = samping_frequency/signal_length*(0:(signal_length/2));

    amplitude = abs(z_frequency/signal_length);
    single_amplitude = amplitude(1:signal_length/2+1,:);
    single_amplitude(2:end-1,:) = 2*single_amplitude(2:end-1,:);
    
    [~,max_index] = max(single_amplitude,[],1);
    fundamental_frequency = frequency(max_index)';

end

function periodicity_error = get_periodicity_error(t,z,fundamental_frequency)
num_signals = size(z,1);

periodicity_error = zeros(num_signals,1);
for iSignal = 1:num_signals
    period = 2*pi/fundamental_frequency(iSignal);
    if period > t(end)
        periodicity_error(iSignal) = inf;
        continue
    end

    z_end = zeros(num_signals,1);
    for jSignal = 1:num_signals
        z_end(jSignal) = interp1(t,z(jSignal,:),period);
    end
    periodicity_error(iSignal) = norm(z_end - z(:,1))/norm(z(:,1));
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