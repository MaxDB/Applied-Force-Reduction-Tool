% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
Parameters.m = 1;
Parameters.k = [1000;1000;30;1;1];
Parameters.L = 1;

%%linear mass
dofs = 3;
eom.M = eye(dofs).*Parameters.m;

%%linear stiffness
k1 = Parameters.k(1);
k2 = Parameters.k(2);
k3 = Parameters.k(3);
k4 = Parameters.k(4);
k5 = Parameters.k(5);
eom.K = [k1 + k2, 0,      0;
         0,      k4 + k5,  0;
         0,         0,   k3];

%%nonlinear restoring force
eom.fnx = @(x) nonlinear_restoring_force(x,Parameters,eom.K);

%%linear damping
% EoM.C

%%applied force
% EoM.F

%%potential energy
eom.V = @(x) potential_energy(x,Parameters);

%% 
system_name = mfilename;
Analytic_Eom = Analytic_System(system_name,eom,"parameters",Parameters);


%-------------------------------------------------------------------------%
function fnx = nonlinear_restoring_force(x,Parameters,K_lin)
k1 = Parameters.k(1);
k2 = Parameters.k(2);
k3 = Parameters.k(3);
k4 = Parameters.k(4);
k5 = Parameters.k(5);

L = Parameters.L;
y2 = x(3);
y1 = x(2);
x1 = x(1);

L_1 = sqrt((x1+L)^2 + y1^2);
L_2 = sqrt((x1-L)^2 + (y1-y2)^2);

L_4 = sqrt(x1^2 + (y1 - L)^2);
L_5 = sqrt(x1^2 + (y1+L)^2);


f1 = k1*(x1+L)*(1-L/L_1) + k2*(x1-L)*(1-L/L_2) + k4*x1*(1-L/L_4) + k5*x1*(1-L/L_5);
f2 = k1*y1*(1-L/L_1) + k2*(y1-y2)*(1-L/L_2) + k4*(y1-L)*(1-L/L_4) + k5*(y1+L)*(1 - L/L_5);
f3 = k2*(y2 - y1)*(1-L/L_2) + k3*y2;

fnx = [f1;f2;f3] - K_lin*x;
% fnx = [f1;f2;f3];
end

function V = potential_energy(x,Parameters)
k1 = Parameters.k(1);
k2 = Parameters.k(2);
k3 = Parameters.k(3);
k4 = Parameters.k(4);
k5 = Parameters.k(5);

L = Parameters.L;
y2 = x(3);
y1 = x(2);
x1 = x(1);

L_1 = sqrt((x1+L)^2 + y1^2);
L_2 = sqrt((x1 - L)^2 + (y1-y2)^2);
L_3 = L + y2;
L_4 = sqrt(x1^2 + (y1 - L)^2);
L_5 = sqrt(x1^2 + (y1+L)^2);

V = 0.5*k1*(L_1 - L)^2 +...
    0.5*k2*(L_2 - L)^2 +...
    0.5*k3*(L_3 - L)^2 +...
    0.5*k4*(L_4 - L)^2 +...
    0.5*k5*(L_5 - L)^2;
end