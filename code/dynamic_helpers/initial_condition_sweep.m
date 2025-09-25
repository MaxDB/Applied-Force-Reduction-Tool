function periodic_solution = initial_condition_sweep(Rom,frequency,ic_limits)
sweep_resolution = 50;
step_reduction = 5.5;
MAX_ERROR = 1e-12;

Model = Rom.Model;
num_modes = size(Model.reduced_modes,2);
period = 2*pi/frequency;
options = odeset(RelTol=1e-8,AbsTol=1e-10);
eom = Rom.get_equation_of_motion();


centre_point = zeros(num_modes,1);

solution_converged = 0;
figure
hold on
while ~solution_converged

    test_conditions = zeros(num_modes,sweep_resolution);
    for iMode = 1:num_modes
        test_conditions(iMode,:) = linspace(centre_point(iMode)-ic_limits(iMode),centre_point(iMode) + ic_limits(iMode),sweep_resolution);
        % test_conditions(iMode,:) = linspace(0,ic_limits(iMode),SWEEP_RESOLUTION);
    end
    num_ics = sweep_resolution^num_modes;


    ic_counter = ones(1,num_modes);
    ic_counter(end) = 0;

    % periodic_sol = [];
    % sol_period = [];
    % sol_periodicity = [];
    % sol_ics = zeros(2*num_modes,0);
    % sol_counter = 0;
    % periodic_solutions = {};

    periodicity_error = zeros(1,num_ics);
    kinetic_energy = zeros(1,num_ics);
    initial_conditions = zeros(num_modes,num_ics);
    max_ke = 0;
    ic_combinations = zeros(num_ics,num_modes);
    for iCondition = 1:num_ics
        ic_counter = increment_counter(ic_counter,sweep_resolution);
        ic_combinations(iCondition,:) = ic_counter;
    end
    
    create_parallel_pool("current");
    test_conditions_const = parallel.pool.Constant(test_conditions);
    parfor iCondition = 1:num_ics
        ic_counter = ic_combinations(iCondition,:);
        initial_condition = zeros(2*num_modes,1);
        for iMode = 1:num_modes
            initial_condition(iMode) = test_conditions_const.Value(iMode,ic_counter(iMode));
        end
        initial_conditions(:,iCondition) = initial_condition(1:num_modes);
        if sum(abs(initial_condition)) == 0
            continue
        end


        Sol = ode45(@(t,z) eom(t,z),[0,period],initial_condition,options);
        periodicity_error(iCondition) = get_periodicity_error(Sol.y(:,end),initial_condition);
        ke = 1/2*sum(Sol.y((num_modes+1):(2*num_modes),:).^2,1);
        max_ke = max(max_ke,max(ke));
        kinetic_energy(iCondition) = ke(end);
    end
    
    % plot3(initial_conditions(1,:),initial_conditions(2,:),periodicity_error,"x");
    plot3(initial_conditions(1,:),initial_conditions(2,:),kinetic_energy,"x");
    ax = gca;
    ax.ZScale = "log";
    drawnow

    % [min_periodicity,min_index] = min(periodicity_error);
    [min_error,min_index] = min(kinetic_energy);
    initial_condition = [initial_conditions(:,min_index);zeros(num_modes,1)];
    
    if min_error < max_ke*MAX_ERROR
        solution_converged = 1;
    end
    ic_step = diff(test_conditions(:,[1,2]),1,2);
    centre_point = initial_condition(1:num_modes);

    if ~any(ismember(abs(centre_point),ic_limits))
        ic_limits = step_reduction*ic_step;
    end
    step_reduction = 2.5;
    sweep_resolution = 11;
end
hold off
num_output_points = 101;
t_sol = linspace(0,period,num_output_points);

[x,y] = ode45(@(t,z) eom(t,z),t_sol,initial_condition,options);
y = y';
x = x';

% 
% y_fourier = zeros(2*num_modes,num_output_points);
% for iMode = 1:2*num_modes
%     y_fourier(iMode,:) = [interpft(y(iMode,1:(end-1)),num_output_points-1),y(iMode,1)];
% end

periodic_solution = {x,y};
% 
% [~,sort_index] = sort(abs(sol_period-period));
% potential_ics = sol_ics(:,sort_index);
% periodic_solutions = periodic_solutions(sort_index);
end

function counter = increment_counter(counter,max_value)
num_counts = size(counter,2);
for iCount = 1:num_counts
    reverse_index = num_counts + 1 -iCount;
    if counter(reverse_index) < max_value
        counter(reverse_index) = counter(reverse_index) + 1;
        break
    else
        counter(reverse_index) = 1;
    end
end

end

%----------
function [periodicity_error,is_terminal,direction] = periodicity_event_fun(t,z,initial_condition,min_period)
is_terminal = 1;
direction = 0;
MAX_ERROR = 1e-4;
%------------------
if t < min_period
    periodicity_error = inf;
    return
end

periodicity_error = get_periodicity_error(z,initial_condition);

if periodicity_error < MAX_ERROR
    periodicity_error = 0;
end
end

%----------
function periodicity_error = get_periodicity_error(zEnd,zStart)

disp_index = 1:(size(zStart,1)/2);
periodicity_error = norm(zEnd(disp_index) - zStart(disp_index))/norm(zStart(disp_index));
end

function [kinetic_energy,is_terminal,direction] = energy_event_fun(t,z,target_period)
is_terminal = 1;
direction = 0;

state_size = size(z,1);
vel = z((state_size/2+1):state_size,:);
kinetic_energy = 1/2* sum(vel.^2);

if t< target_period*3/4
    kinetic_energy = 1;
end

end





% %_---------
% omega_n = sqrt(Rom.Model.reduced_eigenvalues(1));
% period = 2*pi/omega_n;
% max_period = 3*period;
% 
% eom = Rom.get_equation_of_motion();
% input_limits = Rom.Force_Polynomial.input_limit;
% r_lim = [min(input_limits,[],2),max(input_limits,[],2)];
% 
% num_inputs = size(r_lim,1);
% 
% possible_ics = zeros(num_inputs,SWEEP_RESOLUTION);
% for iInput = 1:num_inputs
%     possible_ics(iInput,:) = linspace(r_lim(iInput,1),r_lim(iInput,2),SWEEP_RESOLUTION);
% end
% 
% index = ones(num_inputs,1);
% num_ics = SWEEP_RESOLUTION^num_inputs;
% 
% figure
% hold on
% low_error_ic = [];
% for iIC = 1:num_ics
%     ic = zeros(2*num_inputs,1);
%     for iInput = 1:num_inputs
%         ic(iInput) = possible_ics(iInput,index(iInput));
%     end
% 
%     [t,z] = ode45(@(t,z) eom(t,z),linspace(0,max_period,50),ic);
%     z = z';
%     t = t';
%     fundamental_frequency = check_fundamental_frequency(t,z');
%     periodicity_error = get_periodicity_error(t,z(1:num_inputs,:),fundamental_frequency(1:num_inputs));
%     if all(periodicity_error < MAX_ERROR)
%         low_error_ic(:,end+1) = ic(1:num_inputs)'; %#ok<AGROW>
%     end
%     plot3(ic(1),ic(2),min(periodicity_error),"x")
%     % plot3(ic(1),ic(2),max(periodicity_error),"x")
%     index = increment_index(index,SWEEP_RESOLUTION);
% end
% hold off
% 
% end
% 
% function fundamental_frequency = check_fundamental_frequency(t,z)
%     z_frequency = fft(z);
% 
%     signal_length = size(t,2);
%     samping_frequency = 2*pi/diff(t(1:2));
%     frequency = samping_frequency/signal_length*(0:(signal_length/2));
% 
%     amplitude = abs(z_frequency/signal_length);
%     single_amplitude = amplitude(1:signal_length/2+1,:);
%     single_amplitude(2:end-1,:) = 2*single_amplitude(2:end-1,:);
% 
%     [~,max_index] = max(single_amplitude,[],1);
%     fundamental_frequency = frequency(max_index)';
% 
% end
% 

% 
% 
% function index = increment_index(index,max_index)
% index_size = size(index,1);
% counter = index_size;
% while index(counter) == max_index
%     index(counter) = 1;
%     counter = counter - 1;
%     if counter == 0
%         return
%     end
% end
% index(counter) = index(counter) + 1; 
% end