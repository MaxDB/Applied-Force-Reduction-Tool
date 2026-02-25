% M*x_ddot + C*x_dot + K*x + f(x) = F 
%%Parameters
Parameters.m = 0.1;
Parameters.k1 = 1e3;
Parameters.k3 = 1e3;
% Parameters.c = 1;
Parameters.dofs = 10;

%%linear mass

eom.M = eye(Parameters.dofs).*Parameters.m;

%%linear damping and stiffness
diag_matrix = zeros(Parameters.dofs+2);
diag_index = ((1:Parameters.dofs)'.*[1,1])+1;
diag_right_index = diag_index + [0,1];
diag_left_index = diag_index + [0,-1];

for iDof = 1:Parameters.dofs
    diag_matrix(diag_index(iDof,1),diag_index(iDof,2)) = 2;
    diag_matrix(diag_left_index(iDof,1),diag_left_index(iDof,2)) = -1;
    diag_matrix(diag_right_index(iDof,1),diag_right_index(iDof,2)) = -1;
end
diag_matrix(1,:) = [];
diag_matrix(:,1) = [];
diag_matrix(end,:) = [];
diag_matrix(:,end) = [];

diag_matrix(end,end) = 1; %free BC

eom.K = diag_matrix*Parameters.k1;
% eom.C = diag_matrix*Parameters.c;

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


freq = sqrt(eig(eom.K,eom.M)); 

%-------------------------------------------------------------------------%
function fnx = nonlinear_restoring_force(x,Parameters)
k3 = Parameters.k3;
N = Parameters.dofs;

if class(x) == "double"
fnx = zeros(N,1);
else
fnx = sym("fnx",[N,1]);
end

fnx(1) = x(1)^3 + (x(1)-x(2))^3;
for iDof = 2:(N-1)
    fnx(iDof) = (x(iDof)-x(iDof-1))^3 + (x(iDof)-x(iDof+1))^3;
end
% fnx(N) =  (x(N)-x(N-1))^3 + x(N)^3; %fixed bc
fnx(N) =  (x(N)-x(N-1))^3; %free bc

fnx = k3*fnx;
end

function V = potential_energy(x,Parameters)
k1 = Parameters.k1;
k3 = Parameters.k3;
N = Parameters.dofs;


V = k1/2*x(1)^2 + k3/4*x(1)^4;
for iDof = 1:(N-1)
    V = V + k1/2*(x(iDof)-x(iDof+1))^2 + k3/4*(x(iDof)-x(iDof+1))^4;
end
% V = V+ k1/2*x(N)^2 + k3/4*x(N)^4; %fixed bc

end