function [disp_sep,lambda_sep] = find_sep_rom(rom,force_ratio,target_loadcases)
TARGET_LOADCASES = 100; %approx number of points from origin to end of SEP
% INITIAL_ARC_RADIUS = 1;

MAX_LOADCASES = 1000;    %maximum points per SEP
MAX_INCREMENTS = 1000;   %maximum interations for convergence

CONVERGENCE_TOLERACE = 1e-5;
PSI = 1;

if nargin == 2
    target_loadcases = TARGET_LOADCASES;
end
%-----------------------------%
K = rom.Reduced_Stiffness_Polynomial;
f = rom.Force_Polynomial;
V = rom.Potential_Polynomial;

num_modes = size(f,1);



%start with linear solution
r_0 = zeros(num_modes,1);
lambda_0 = 0;

target_lambda_inc = 1/target_loadcases;
base_force = force_ratio;
% arc_radius = INITIAL_ARC_RADIUS;

%initial arc radius
lin_stiffness = K.evaluate_polynomial(r_0);
r_lin = r_0 + lin_stiffness\(target_lambda_inc*base_force);
lin_force = lin_stiffness*r_lin;
lambda_lin = (base_force/target_loadcases) ./ lin_force;
lambda_lin = mean(lambda_lin(~isnan(lambda_lin)));
arc_radius = (r_lin'*r_lin + PSI^2*lambda_lin^2*(lin_force'*lin_force)/target_loadcases^2);

lambda_sep = zeros(1,MAX_LOADCASES);
disp_sep = zeros(num_modes,MAX_LOADCASES);

for iLoad = 1:MAX_LOADCASES
    lambda = lambda_0 + target_lambda_inc;

    lin_stiffness = K.evaluate_polynomial(r_0);
    r = r_0 + lin_stiffness\(target_lambda_inc*base_force);

    

    control_matrix = zeros(num_modes+1);
    mode_span = 1:num_modes;
    for iInc = 1:MAX_INCREMENTS
        % Kt = K.evaluate_polynomial(r_0+r);
        % f_x = f.evaluate_polynomial(r_0+r);
        Kt = K.evaluate_polynomial(r);
        f_x = f.evaluate_polynomial(r);
        % eq_condition = -f_x + (lambda_0+lambda)*base_force;
         eq_condition = -f_x + lambda*base_force;
        

        
        % convergence_test = eq_condition./((lambda_0+lambda)*base_force.*f_x);
        convergence_test = eq_condition./((lambda*base_force).*f_x);
        if max(abs(convergence_test(~isinf(convergence_test)))) < CONVERGENCE_TOLERACE
            break
        end

        % arc_corrector =  arc_radius^2 - (r'*r + PSI^2*lambda^2*(base_force'*base_force));
         arc_corrector =  arc_radius^2 - ((r-r_0)'*(r-r_0) + PSI^2*(lambda-lambda_0)^2*(base_force'*base_force));

        control_matrix(mode_span,mode_span) = Kt;
        control_matrix(mode_span,num_modes+1) = -base_force;
        control_matrix(num_modes+1,mode_span) = 2*(r-r_0)';
        control_matrix(num_modes+1,num_modes+1) = 2*PSI^2*(lambda-lambda_0)*(base_force'*base_force);

        increment = control_matrix\[eq_condition;arc_corrector];

        r = r + increment(1:num_modes);
        lambda = lambda + increment(end);
    end
    if iInc == MAX_INCREMENTS
        warning("Convergence tolerance not met for arc length continuation of SEPs")
        break
    end
    
   
    % r_0 = r_0 + r;
    % lambda_0 = lambda_0 + lambda;
    r_0 = r;
    lambda_diff = lambda-lambda_0;
    lambda_0 = lambda;

    if V.evaluate_polynomial(r_0) > rom.Model.fitting_energy_limit
        break
    end

    lambda_sep(1,iLoad) = lambda_0;
    disp_sep(:,iLoad) = r_0;

    
    lambda_ratio = lambda_diff/target_lambda_inc;
    if lambda_ratio > 1.1
        arc_radius = arc_radius/min(lambda_ratio,2);
    elseif lambda_ratio < 1/1.1
        arc_radius = arc_radius*min(1/lambda_ratio,2);
    end
end
lambda_sep(:,iLoad:end) = [];
disp_sep(:,iLoad:end) = [];
end