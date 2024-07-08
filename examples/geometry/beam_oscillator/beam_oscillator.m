% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
% omega_nr = pi^2;
% omega_ns = 4*pi^2;
% psi_3_1 = 0.35;
% psi_1_3 = 0.35;
% psi_2_2 = 2*pi^4;
% psi_4_0 = pi^4/2;
% psi_0_4 = 8*pi^4;

omega_nr = 1;
omega_ns = 4;
psi_3_1 = 0.35; %a2
psi_1_3 = 0;    %a4
psi_2_2 = 2;    %a3
psi_4_0 = 1/2;  %a1
psi_0_4 = 8;    %a5

%%linear mass
eom.M = eye(2);

%%linear stiffness
eom.K = [omega_nr^2,0;0,omega_ns^2];

%%nonlinear restoring force
eom.f = @(x) nonlinear_restoring_force(x,psi_3_1,psi_1_3,psi_2_2,psi_4_0,psi_0_4);

%%linear damping
% EoM.C

%%applied force
% EoM.F

%%potential energy
eom.V = @(x) potential_energy(x,omega_nr,omega_ns,psi_3_1,psi_1_3,psi_2_2,psi_4_0,psi_0_4);

%% 
system_name = mfilename;
Analytic_Eom = Analytic_System(system_name,eom);


%-------------------------------------------------------------------------%
function fnx = nonlinear_restoring_force(x,psi_3_1,psi_1_3,psi_2_2,psi_4_0,psi_0_4)
f_r = psi_4_0*x(1).^3 + 3*psi_3_1*x(1).^2.*x(2) + psi_2_2*x(1).*x(2).^2 + psi_1_3*x(2).^3;
f_s = psi_3_1*x(1).^3 + psi_2_2*x(1).^2.*x(2) + 3*psi_1_3*x(1).*x(2).^2 + psi_0_4*x(2).^3;

fnx = [f_r;f_s];
end

function V = potential_energy(x,omega_nr,omega_ns,psi_3_1,psi_1_3,psi_2_2,psi_4_0,psi_0_4)
V = 0.5*omega_nr^2*x(1).^2 + 0.5*omega_ns^2*x(2).^2 + 0.25*psi_4_0*x(1).^4  ...
+ psi_3_1*x(1).^3.*x(2) + 0.5*psi_2_2*x(1).^2.*x(2).^2 + psi_1_3*x(1).*x(2).^3 + 0.25*psi_0_4*x(2).^4;
end