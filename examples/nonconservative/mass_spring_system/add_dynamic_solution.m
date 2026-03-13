clear
% close all
set_visualisation_level(1)

system_name = "mass_spring_roller_0";
Dyn_Data = initalise_dynamic_data(system_name);

Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 2;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e-1;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;

Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 8;


%-----------------------------------------%
Continuation_Opts.energy_limit_multiplier = 1.2;
% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts,"type","fom");
% Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts,"type","rom");

%-----------------------------------
Model = Dyn_Data.Dynamic_Model.Model;
[mass_factor,stiffness_factor] = get_damping_coeffs(Model,1,0.005);
Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = mass_factor;
Damping_Data.stiffness_factor = stiffness_factor;


Force_Data.type = "shape";
Force_Data.shape = [0;1;0];
Force_Data.continuation_variable = "frequency";
Force_Data.frequency = 1.1;
Force_Data.amplitude = 7e-3;


% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e1;
Continuation_Opts.min_inc = 1e-3;

Continuation_Opts.forward_steps = 10;
Continuation_Opts.backward_steps = 10;
Continuation_Opts.frequency_range = [0.995,3.2];
%-----------------------------------------%
% 

% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","rom");
% Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","fom");


orbit = Dyn_Data.get_orbit(1,133);
Model = Dyn_Data.Dynamic_Model.Model;
Force_Data.frequency = 2*pi/orbit.T;
eom = Model.get_equation_of_motion("damping",Damping_Data,"forcing",Force_Data);


t0 = 0;
x0 = orbit.xbp(1,:)';

[t,y] = ode45(eom,[0,5000],x0);
y=y';
period = 2*pi/Force_Data.frequency;
num_periods = t(end)/period;
period_range = floor(num_periods) - [1,0];
time_range = period_range*period;
orbit_index = t > time_range(1) & t < time_range(2);
last_index = find(orbit_index,1,"last");
orbit_index(last_index+1) = true;
y_orbit = y(:,orbit_index);
t_orbit = t(orbit_index);
t_orbit = t_orbit - period*period_range(1);
initial_conditions = {t_orbit,y_orbit,Force_Data.frequency};

Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 1e1;
Continuation_Opts.min_inc = 1e-3;

Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 170;
% Continuation_Opts.collocation_degree = 12;
% Continuation_Opts.min_discretisation_num = 40;
Continuation_Opts.energy_limit_multiplier = 6;

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts,"type","fom","ic",initial_conditions);


function [alpha,beta] = get_damping_coeffs(Model,modes,damping_ratios)
[evec,eval] = eig(Model.stiffness,Model.mass);
num_modes = length(modes);

freq = sqrt(diag(eval));

if num_modes == 1
    alpha = 0;
    beta = 2*damping_ratios/freq(modes);
else

    coeff_matrix = [ones(num_modes,1),freq(modes).^2];
    coeffs = 2*damping_ratios*coeff_matrix\freq(modes);

    alpha = coeffs(1);
    beta = coeffs(2);
end
end



%moderate damping
plot_index = 1:5;
compare_solutions("energy","mass_spring_roller_0",plot_index,"mass_spring_roller_1",plot_index,"mass_spring_roller_12",plot_index);


% %moderate damping validation
% plot_index = 2;
% compare_solutions("energy","mass_spring_roller_1",plot_index,"mass_spring_roller_12",plot_index,"validation",[1,0]);

%small damping
plot_index = [1,6:11];
compare_solutions("energy","mass_spring_roller_0",plot_index,"mass_spring_roller_12",plot_index);

compare_solutions("energy","mass_spring_roller_0",plot_index,"mass_spring_roller_1",plot_index);

%small damping 7e-3 + 8e-3 amp
plot_index = [1,10:12];
compare_solutions("energy","mass_spring_roller_0",plot_index);

compare_solutions("energy","mass_spring_roller_0",plot_index,"mass_spring_roller_1",plot_index);

%small damping validation
compare_solutions("energy","mass_spring_roller_0",[10,12],"mass_spring_roller_1",10,"validation",[0,1]);


plot_index = [1,5];
compare_solutions("energy","mass_spring_roller_0",plot_index,"mass_spring_roller_1",plot_index,"mass_spring_roller_12",plot_index);

compare_validation("mass_spring_roller_1","energy",5,2);