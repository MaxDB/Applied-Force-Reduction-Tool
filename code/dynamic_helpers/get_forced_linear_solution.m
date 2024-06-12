function [t0,z0,mode_amplitude] =  get_forced_linear_solution(Rom,Nonconservative_Input,type)
TIME_RESOLUTION = 100; %number of time points in solution
INITIAL_FORCE = 0.01; %initual guess at force amplitude for a linear response 
MAX_ITERATIONS = 100; %maximum number of attempts to find linear solution
MAX_LINEARITY_ERROR = 1; %maximum allowable percentage difference in approximate linear solution

%compare to restoring force to get system independent guess/

Model = Rom.Model;
%find linear solution
r_eigenvectors = Model.reduced_eigenvectors;
r_eigenvalues = Model.reduced_eigenvalues;
r_modes = Model.reduced_modes;
num_modes = size(r_modes,2);

damping = r_eigenvectors'*Nonconservative_Input.damping*r_eigenvectors;


mode_map = Nonconservative_Input.mode_map;
mode_index = find(mode_map);
if isempty(mode_map)
    error("Mode: " + backbone_num + " not captured by ROM. Add the corresponding mode.")
end


mode_damping = damping(mode_map,mode_map);

frequency = Nonconservative_Input.frequency;
period = 2*pi/frequency;
mode_natural_frequency = sqrt(r_eigenvalues(mode_map,mode_map));

freq_squared_diff = mode_natural_frequency^2 - frequency^2;
amp_denominator = freq_squared_diff^2 + mode_damping^2*frequency^2;
sin_amp_scale_factor = freq_squared_diff/amp_denominator;
cos_amp_scale_factor = -mode_damping*frequency/amp_denominator;

t = linspace(0,period,TIME_RESOLUTION);
t_prod = frequency*t;
sin_t = sin(t_prod);
cos_t = cos(t_prod);

Eom_Input = Rom.get_solver_inputs("coco_frf",Nonconservative_Input);

input_order = Eom_Input.input_order;
Force_Data = Eom_Input.Force_Data;
Disp_Data = Eom_Input.Disp_Data;
Damping_Data = Eom_Input.Damping_Data;
Applied_Force_Data = Eom_Input.Applied_Force_Data;

eom = @(z,amp) coco_forced_eom(t,z,amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data);



mode_amplitude = INITIAL_FORCE;
for iIteration = 1:MAX_ITERATIONS
    sin_amp = mode_amplitude*sin_amp_scale_factor;
    cos_amp = -mode_amplitude*cos_amp_scale_factor;

    r_lin = sin_amp*sin_t + cos_amp*cos_t;
    r_dot_lin = frequency*sin_amp*cos_t - frequency*cos_amp*sin_t;
    r_ddot_lin = -(frequency^2)*sin_amp*sin_t - (frequency^2)*cos_amp*cos_t;
    
    x_lin = zeros(2*num_modes,TIME_RESOLUTION);
    x_lin(mode_index,:) = r_lin;
    x_lin(mode_index + num_modes,:) = r_dot_lin;
    z_dot_sol = eom(x_lin,mode_amplitude);
    
    r_ddot_sol = z_dot_sol(mode_index+num_modes,:);
    r_ddot_sol(abs(r_ddot_sol) < 1e-10) = 0;
    acceleration_diff = abs(r_ddot_lin - r_ddot_sol)./max(abs(r_ddot_lin));
    max_error = max(abs(acceleration_diff(~isinf(acceleration_diff))));
    if max_error < MAX_LINEARITY_ERROR/100
        break
    elseif iIteration == MAX_ITERATIONS
        error("Could not converge to linear solution")
    else
        mode_amplitude = mode_amplitude/2;
    end
end
z0 = x_lin;
t0 = t;
end