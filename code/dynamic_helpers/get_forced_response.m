function [t_sol,z_sol] = get_forced_response(Rom,Nonconservative_Input,period)
MAX_INCREMENTS = 100;
MAX_ERROR = 1e-4;

%-------------------------------------------------------------------------%
Eom_Input = Rom.get_solver_inputs("coco_frf",Nonconservative_Input);

input_order = Eom_Input.input_order;
Force_Data = Eom_Input.Force_Data;
Disp_Data = Eom_Input.Disp_Data;
Damping_Data = Eom_Input.Damping_Data;
Applied_Force_Data = Eom_Input.Applied_Force_Data;

amp = Applied_Force_Data.amplitude; 


eom = @(t,z) coco_forced_eom(t,z,amp,period,input_order,Force_Data,Disp_Data,Damping_Data,Applied_Force_Data);
%-------------------------------------------------------------------------%
num_r_modes = size(Rom.Model.reduced_modes,2);


z0 = zeros(2*num_r_modes,1);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
t = 0;
t_increment = period;

for iInc = 1:MAX_INCREMENTS
    t_span = t(end) + [0,t_increment];
    [t,z] = ode45(eom,t_span,z0,opts);
    z0 = z(end,:);

    error = check_periodicity(z);
    if error < MAX_ERROR
        break
    end
end

if iInc == MAX_INCREMENTS
    error("Could not converge to forced solution")
end

t_sol = t - t(1);
z_sol = z';

function error = check_periodicity(z)
    z_0 = z(1,:);
    z_T = z(end,:);
    error = norm(z_T - z_0)./norm(z_0);
end

end

