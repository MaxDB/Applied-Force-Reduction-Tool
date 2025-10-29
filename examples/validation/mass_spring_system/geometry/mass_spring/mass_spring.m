% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
Parameters.m = 1;
Parameters.k = [0.25;0.75;0.75;1.5;250;750]; %kx1 kx2 ky1 ky2 kz1 kz2
Parameters.L = 1;

%%linear mass
dofs = 3;
eom.M = eye(dofs).*Parameters.m;

%%linear stiffness
kx1 = Parameters.k(1);
kx2 = Parameters.k(2);
ky1 = Parameters.k(3);
ky2 = Parameters.k(4);
kz1 = Parameters.k(5);
kz2 = Parameters.k(6);
L = Parameters.L;
eom.K = [kx1 + kx2, 0,      0;
         0,      ky1 + ky2,  0;
         0,         0,   kz1 + kz2];

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
kx1 = Parameters.k(1);
kx2 = Parameters.k(2);
ky1 = Parameters.k(3);
ky2 = Parameters.k(4);
kz1 = Parameters.k(5);
kz2 = Parameters.k(6);

L = Parameters.L;
z = x(3);
y = x(2);
x = x(1);

L_x1 = sqrt((x-L)^2 + y^2 + z^2);
L_y1 = sqrt(x^2 + (y-L)^2 + z^2);
L_z1 = sqrt(x^2 + y^2 + (z-L)^2);

L_x2 = sqrt((x+L)^2 + y^2 + z^2);
L_y2 = sqrt(x^2 + (y+L)^2 + z^2);
L_z2 = sqrt(x^2 + y^2 + (z+L)^2);

f1 = kx1*(x-L)*(1-L/L_x1) + ky1*x*(1-L/L_y1) + kz1*x*(1-L/L_z1) + kx2*(x+L)*(1-L/L_x2) + ky2*x*(1-L/L_y2) + kz2*x*(1-L/L_z2);
f2 = kx1*y*(1-L/L_x1) + ky1*(y-L)*(1-L/L_y1) + kz1*y*(1-L/L_z1) + kx2*y*(1-L/L_x2) + ky2*(y+L)*(1-L/L_y2) + kz2*y*(1-L/L_z2);
f3 = kx1*z*(1-L/L_x1) + ky1*z*(1-L/L_y1) + kz1*(z-L)*(1-L/L_z1) + kx2*z*(1-L/L_x2) + ky2*z*(1-L/L_y2) + kz2*(z+L)*(1-L/L_z2);

fnx = [f1;f2;f3] - K_lin*[x;y;z];
% fnx = [f1;f2;f3];
end

function V = potential_energy(x,Parameters)
kx1 = Parameters.k(1);
kx2 = Parameters.k(2);
ky1 = Parameters.k(3);
ky2 = Parameters.k(4);
kz1 = Parameters.k(5);
kz2 = Parameters.k(6);
L = Parameters.L;
z = x(3);
y = x(2);
x = x(1);

L_x1 = sqrt((x-L)^2 + y^2 + z^2);
L_y1 = sqrt(x^2 + (y-L)^2 + z^2);
L_z1 = sqrt(x^2 + y^2 + (z-L)^2);

L_x2 = sqrt((x+L)^2 + y^2 + z^2);
L_y2 = sqrt(x^2 + (y+L)^2 + z^2);
L_z2 = sqrt(x^2 + y^2 + (z+L)^2);

V = 0.5*kx1*(L_x1 - L)^2 + 0.5*kx2*(L_x2 - L)^2 +...
    0.5*ky1*(L_y1 - L)^2 + 0.5*ky2*(L_y2 - L)^2 +...
    0.5*kz1*(L_z1 - L)^2 + 0.5*kz2*(L_z2 - L)^2;
end