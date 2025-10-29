% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
Parameters.omega(1) = 1;
Parameters.omega(2) = 5;
%%linear mass
dofs = 2;
Eom.M = eye(dofs);

%%linear stiffness
omega = Parameters.omega;
Eom.K = diag(omega.^2);

%%nonlinear restoring force
Eom.fnx = @(x) nonlinear_restoring_force(x,Parameters);

%%linear damping
% EoM.C

%%applied force
% EoM.F

%%potential energy
Eom.V = @(x) potential_energy(x,Parameters);

%% 
system_name = mfilename;
Analytic_Eom = Analytic_System(system_name,Eom,"parameters",Parameters);


%-------------------------------------------------------------------------%
function fnx = nonlinear_restoring_force(x,Parameters)
omega_1 = Parameters.omega(1);
omega_2 = Parameters.omega(2);

x1 = x(1);
x2 = x(2);

f1  = omega_1^2/2*(3*x1^2 + x2^2) + omega_2^2*x1*x2 + (omega_1^2 + omega_2^2)/2*x1*(x1^2+x2^2);
f2  = omega_2^2/2*(3*x2^2 + x1^2) + omega_1^2*x1*x2 + (omega_1^2 + omega_2^2)/2*x2*(x1^2+x2^2);
fnx = [f1;f2];
end

function V = potential_energy(x,Parameters)
omega_1 = Parameters.omega(1);
omega_2 = Parameters.omega(2);
a1 = omega_1^2/2;
a2 = omega_2^2/2;

x1 = x(1);
x2 = x(2);

V = a1*x1^2 + a2*x2^2 + a1*x1^3 + a2*x2^3 + ...
    a1*x2^2*x1 + a2*x1^2*x2 + ...
    (a1+a2)*(x1^4/4 + x2^4/4 + x1^2*x2^2/2);

end