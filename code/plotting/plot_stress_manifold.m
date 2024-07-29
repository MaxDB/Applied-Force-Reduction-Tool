function plot_stress_manifold(Rom,L_modes)
PLOT_RESOLUTION = 101;


data_path = Rom.data_path;

load_static_data_start = tic;
load(data_path + "Static_Data.mat","Static_Data")

load_static_data_time = toc(load_static_data_start);
log_message = sprintf("Static dataset loaded: %.1f seconds" ,load_static_data_time);
logger(log_message,3)

Rom = Reduced_System(Static_Data,Static_Data.validated_degree,"full");

if L_modes == 0
    Rom.Condensed_Displacement_Polynomial.plot_polynomial
    return
end



Static_Data = Static_Data.add_validation_data(L_modes);
Rom = Reduced_System(Static_Data,Static_Data.validated_degree,"full");
num_dofs = Rom.Model.num_dof;

r_limits = Rom.Condensed_Displacement_Polynomial.input_limit;
r = linspace(r_limits(1),r_limits(2),PLOT_RESOLUTION);
h = linspace(r_limits(1),r_limits(2),PLOT_RESOLUTION);

h_r = h;
h_L = h;

plot_r = zeros(PLOT_RESOLUTION);
plot_h = zeros(PLOT_RESOLUTION);
plot_x = zeros(PLOT_RESOLUTION,PLOT_RESOLUTION,num_dofs);
plot_energy = zeros(PLOT_RESOLUTION);

for iR = 1:PLOT_RESOLUTION
    for iH = 1:PLOT_RESOLUTION
        plot_r(iR,iH) = r(1,iR) + h_r(1,iH);
        plot_h(iR,iH) = h_L(1,iH);
        h_i = [h_r(1,iH);h_L(1,iH)];
        plot_x(iR,iH,:) = Rom.expand(r(1,iR),h_i);
        r_energy = Rom.Potential_Polynomial.evaluate_polynomial(r(1,iR));
        h_energy = 1/2*((Rom.Low_Frequency_Stiffness_Polynomial.evaluate_polynomial(r(1,iR)))*h_i)'*h_i;
        plot_energy(iR,iH) = r_energy + h_energy;
    end
end


energy_limit = Rom.Model.energy_limit;
remove_index = plot_energy>energy_limit;
plot_r(remove_index) = nan;
plot_h(remove_index) = nan;

h_1d = zeros(1,PLOT_RESOLUTION);
x_1d = Rom.expand(r);

for iDof = 1:num_dofs
    
    figure
    mesh(plot_r,plot_h,plot_x(:,:,iDof))
    box on
    xlabel("r + h_r")
    ylabel("h_L")
    zlabel("x_{" + iDof + "}")
    
    hold on
    plot3(r,h_1d,x_1d(iDof,:),'k')
    hold off
end






end