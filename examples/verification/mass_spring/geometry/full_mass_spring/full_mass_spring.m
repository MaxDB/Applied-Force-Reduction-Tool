% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
Parameters.omega(1) = 1;
Parameters.omega(2) = 5;
%%linear mass
dofs = 2;
Eom.M = eye(dofs);

%%linear stiffness
omega = Parameters.omega;
% Eom.K = diag(omega.^2);

%%nonlinear restoring force
Eom.f = @(x) restoring_force(x,Parameters);

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
function f = restoring_force(x,Parameters)
omega_1 = Parameters.omega(1);
omega_2 = Parameters.omega(2);

x1 = x(1);
x2 = x(2);

lambda1 = (1 - 1/sqrt(x2^2 + (1+x1)^2));
lambda2 = (1 - 1/sqrt(x1^2 + (1+x2)^2));

f1 = omega_1^2*(x1+1)*lambda1 + omega_2^2*x1*lambda2;
f2 = omega_1^2*x2*lambda1 + omega_2^2*(x2+1)*lambda2;

f = [f1;f2];
end

function V = potential_energy(x,Parameters)
omega_1 = Parameters.omega(1);
omega_2 = Parameters.omega(2);

x1 = x(1);
x2 = x(2);

L1 = sqrt(x2^2+(1+x1)^2);
L2 = sqrt(x1^2+(1+x2)^2);
V1 = 0.5*omega_1^2*(L1-1)^2;
V2 = 0.5*omega_2^2*(L2-1)^2;

V = V1+V2;
end