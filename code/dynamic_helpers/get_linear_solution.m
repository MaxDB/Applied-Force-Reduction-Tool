function [t0,z0] = get_linear_solution(Rom,backbone_num,type)
TIME_RESOLUTION = 100; %number of time points in solution
INITIAL_AMP_SCALE_FACTOR = 0.01; %initial guess at ratio of end of linear regieme to max displacement
MAX_ITERATIONS = 100; %maximum number of attempts to find linear solution
MAX_LINEARITY_ERROR = 1e-3; %maximum allowable percentage difference in approximate linear solution

Model = Rom.Model;
%find linear solution
r_eigenvalues = Model.reduced_eigenvalues;
r_modes = Model.reduced_modes;


mode_index = find(r_modes == backbone_num);


natural_freq = sqrt(r_eigenvalues(mode_index));
period = 2*pi/natural_freq;
t = linspace(0,period,TIME_RESOLUTION);
zeta = zeros(1,TIME_RESOLUTION);

switch type
    case "rom"
        if isempty(mode_index)
            error("Backbone " + backbone_num + " not captured by ROM. Add the corresponding mode.")
        end

        num_modes = length(r_modes);
        y = zeros(2*num_modes,TIME_RESOLUTION);
        y(mode_index,:) = sin(natural_freq*t);
        y(mode_index+num_modes,:) = natural_freq*cos(natural_freq*t);
        x_ddot = -natural_freq^2*sin(natural_freq*t);
        x_ddot(abs(x_ddot) < 1e-10) = 0;


        Eom_Input = Rom.get_solver_inputs("coco_backbone");

        input_order = Eom_Input.input_order;
        Force_Data = Eom_Input.Force_Data;
        Disp_Data = Eom_Input.Disp_Data;

        eom = @(z) coco_eom(t,z,zeta,input_order,Force_Data,Disp_Data);
    case "fom"
        Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
        num_modes = Model.num_dof;
        y = zeros(2*num_modes,TIME_RESOLUTION);
        y(backbone_num,:) = sin(natural_freq*t);
        y(backbone_num+num_modes,:) = natural_freq*cos(natural_freq*t);
        x_ddot = -natural_freq^2*sin(natural_freq*t);
        x_ddot(abs(x_ddot) < 1e-10) = 0;
        
        Eom_Input = Analytic_Eom.get_solver_inputs("free");
        eom = @(z) direct_eom(0,z,zeta,Eom_Input.modal_restoring_force);
end


amplitude = min(abs(Rom.reduced_displacement_limits(mode_index,:)))*INITIAL_AMP_SCALE_FACTOR;
for iIteration = 1:MAX_ITERATIONS
    z = y*amplitude;
    x_ddot_linear = x_ddot*amplitude;

    z_dot = eom(z);
    x_ddot_sol = z_dot(mode_index+num_modes,:);
    x_ddot_sol(abs(x_ddot_sol) < 1e-10) = 0;
    acceleration_diff = abs(x_ddot_linear - x_ddot_sol)./max(abs(x_ddot_linear));
    max_error = max(abs(acceleration_diff(~isinf(acceleration_diff))));
    if max_error < MAX_LINEARITY_ERROR/100
        break
    elseif iIteration == MAX_ITERATIONS
        error("Could not converge to linear solution")
    else
        amplitude = amplitude/2;
    end
end
z0 = z;
t0 = t;
end