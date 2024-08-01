clear
% M*x_ddot + C*x_dot + K*x + f(x) = F 
system_name = "h_oscillator";
%%Parameters
m = 1;
k1 = 1000;
k2 = 10;
k3 = 1000;
k4 = 1;
k5 = 21.2;
L0 = 1;

%%linear mass
eom.M = eye(4);

%%linear stiffness
eom.K = [k4,0,0,0;
         0,k5,0,0;
         0,0,k1+k2,-k2;
         0,0,-k2,k2+k3];

%%nonlinear restoring force
eom.f = @(x) nonlinear_restoring_force(x,k1,k2,k3,k4,k5,L0);

%%linear damping
% EoM.C

%%applied force
% EoM.F

%%potential energy
eom.V = @(x) potential_energy(x,k1,k2,k3,k4,k5,L0);

%% 
analytic_EoM = Analytic_System(system_name,eom,"save","geometry\" + system_name);


%-------------------------------------------------------------------------%
function fnx = nonlinear_restoring_force(x,k1,k2,k3,k4,k5,L0)
d1 = sqrt(x(1).^2+(L0+x(3)).^2);
d2 = sqrt((x(1)-x(2)).^2+(L0-x(3)+x(4)).^2);
d3 = sqrt(x(2).^2+(L0-x(4)).^2);
d4 = sqrt((L0+x(1)).^2+x(3).^2);
d5 = sqrt((L0+x(2)).^2+x(4).^2);

fy1 = -(k4*L0*(L0+x(1)))/d4 - (k1*x(1)*(L0-d1))/d1 - (k2*(x(1)-x(2))*(L0-d2))/d2 + k4*L0;
fy2 = -(k5*L0*(L0+x(2)))/d5 - (k3*x(2)*(L0-d3))/d3 - (k2*(x(2)-x(1))*(L0-d2))/d2 + k5*L0;
fx1 = -(k4*x(3)*(L0-d4))/d4 - (k1*L0*(L0+x(3)))/d1 - (k2*L0*(x(3)-x(4)-L0))/d2 + (k1-k2)*L0;
fx2 = -(k5*x(4)*(L0-d5))/d5 - (k3*L0*(x(4)-L0))/d3 - (k2*L0*(L0+x(4)-x(3)))/d2 + (k2-k3)*L0;

fnx = [fy1;fy2;fx1;fx2];
end

function V = potential_energy(x,k1,k2,k3,k4,k5,L0)
dL1 = sqrt((L0 + x(3)).^2 + x(1).^2) - L0;
dL2 = sqrt((L0 + x(4) - x(3)).^2 + (x(2)-x(1)).^2) - L0;
dL3 = sqrt((L0 - x(4)).^2 + x(2).^2) - L0;
dL4 = sqrt(x(3).^2 + (L0 + x(1)).^2) - L0;
dL5 = sqrt(x(4).^2 + (L0 + x(2)).^2) - L0;

V = 0.5*(k1*dL1.^2 + k2*dL2.^2 + k3*dL3.^2 + k4*dL4.^2 + k5*dL5.^2);
end