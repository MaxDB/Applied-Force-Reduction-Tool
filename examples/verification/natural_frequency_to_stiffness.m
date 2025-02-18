clear
close all

x_spring_stiffness = 1e3;
target_natural_frequency = [138.287,142.287,146.287,150.287,154.287,158.287,162.287];

min_y_spring_stiffness = 1e3;
max_y_spring_stiffness = 1e5;
num_points = 20;


%--------- Software Settings ---------%
set_logging_level(1)
set_visualisation_level(0)
%-------------------------------------%
Static_Opts = struct([]);
Calibration_Opts.disable_calibration = 1;

%--------- System Settings ---------%
system_name = "ic_demo";
initial_modes = [1];
energy_limit = 2.5;

y_spring_stiffness = logspace(log10(min_y_spring_stiffness),log10(max_y_spring_stiffness),num_points);
natural_frequency = zeros(1,num_points);
for iK = 1:num_points
    set_bc_stiffness(system_name,x_spring_stiffness,y_spring_stiffness(iK))
    Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);
    natural_frequency(iK) = sqrt(Model.reduced_eigenvalues);
end

figure
semilogx(y_spring_stiffness,natural_frequency,"x-")

target_stiffness = interp1(natural_frequency,y_spring_stiffness,target_natural_frequency);