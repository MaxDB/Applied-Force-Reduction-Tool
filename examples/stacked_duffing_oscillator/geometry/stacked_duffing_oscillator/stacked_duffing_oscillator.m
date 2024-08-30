% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
Parameters.m = [0.1;10;10];
% Parameters.k = [2;1;2;1];
% Parameters.alpha = [0;0;1;0];
Parameters.k = [1,2,2,1];
Parameters.alpha = [1;0;0;0];


%%linear mass
dofs = size(Parameters.m,1);
eom.M = eye(dofs).*Parameters.m;

%%linear stiffness
k1 = Parameters.k(1);
k2 = Parameters.k(2);
k3 = Parameters.k(3);
k4 = Parameters.k(4);
eom.K = [k1+k2, -k2,      0;
         -k2,      k2+k3,  -k3;
         0,         -k3,   k3 + k4];

%%nonlinear restoring force
eom.fnx = @(x) nonlinear_restoring_force(x,Parameters);

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
function fnx = nonlinear_restoring_force(x,Parameters)
alpha = Parameters.alpha;
% f1 = alpha(1)*x(1,:).^3 - alpha(2)*(x(2,:) - x(1,:)).^3;
% f2 = alpha(2)*(x(2,:) - x(1,:)).^3 - alpha(3)*(x(3,:) - x(2,:)).^3;
% f3 = alpha(3)*(x(3,:) - x(2,:)).^3;
f1 = alpha(1)*x(1,:).^3 - alpha(2)*(x(2,:) - x(1,:)).^3;
f2 = alpha(2)*(x(2,:) - x(1,:)).^3 - alpha(3)*(x(3,:) - x(2,:)).^3;
f3 = alpha(3)*(x(3,:) - x(2,:)).^3 + alpha(4)*x(3,:).^3;

fnx = [f1;f2;f3];
end

function V = potential_energy(x,Parameters)
alpha = Parameters.alpha;
k = Parameters.k;
V = 0.5*k(1)*x(1).^2 + 0.25*alpha(1)*x(1).^4 + ...
    0.5*k(2)*(x(2)-x(1)).^2 + 0.25*alpha(2)*(x(2)-x(1)).^4 + ...
    0.5*k(3)*(x(3)-x(2)).^2 + 0.25*alpha(3)*(x(3)-x(2)).^4 + ...
    0.5*k(4)*x(3).^2 + 0.25*alpha(4)*(x(3)).^4;
end