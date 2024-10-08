function [t_sol,z_sol] = get_forced_response(Rom,Nonconservative_Input,period)
MAX_INCREMENTS = 100000;
MAX_ERROR = 1e-4;
CONVERGENCE_SPAN = 10;
NUM_PERIODS = 1;

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

if isfield(Nonconservative_Input,"z0")
    z0 = Nonconservative_Input.z0*1.01;
else
    z0 = zeros(2*num_r_modes,1);
end
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
t = 0;
t_increment = NUM_PERIODS*period;
% t_all = [];
% z_all = [];
figure
box on 
hold on
ax = gca;
ax.YScale = "log";
error = zeros(1,CONVERGENCE_SPAN);
mean_error = nan;
for iInc = 1:MAX_INCREMENTS
    t_span = t(end) + [0,t_increment];
    [t,z] = ode45(eom,t_span,z0,opts);
    % t_all = [t_all,t'];
    % z_all = [z_all,z'];
    
    z0 = z(end,:);
    error_index = mod(iInc-1,CONVERGENCE_SPAN)+1;
    error(error_index) = check_periodicity(z(:,1:num_r_modes));
    semilogy(iInc,error(error_index),'.')

    if iInc >= CONVERGENCE_SPAN
        new_mean_error = mean(error);
        semilogy([iInc-1,iInc],[mean_error,new_mean_error],"k","LineWidth",2)
        error_diff = 2*abs(mean_error - new_mean_error)/(mean_error + new_mean_error);
        mean_error = new_mean_error;
    end

    

    drawnow
    if mean_error < MAX_ERROR
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

