clear
close all
load("data\exhaust_17\Static_Data.mat");
%-------------------------------------------------------------------------%
Static_Data = Static_Data.add_validation_data(5);

Rom = Reduced_System(Static_Data);
Disp_Poly = Rom.Physical_Displacement_Polynomial;
H_Displacement_Poly = Rom.Low_Frequency_Displacement_Polynomial;
Force_Poly = Rom.Force_Polynomial;
H_Force_Poly = Rom.Low_Frequency_Force_Polynomial;


%static data
NUM_OUTPUTS = 4;
outputs = randi(size(Disp_Poly,1),1,NUM_OUTPUTS);

r = Static_Data.reduced_displacement;
h = zeros(size(r));
x = Static_Data.physical_displacement(outputs,:);
f = Static_Data.restoring_force;
perturbation = squeeze(Static_Data.condensed_perturbation(:,2,:));

h_displacement = Static_Data.get_h_displacement;
h_L_displacement = squeeze(h_displacement(2,2,:))';

f_L = Static_Data.perturbation_scale_factor*ones(1,size(f,2));


%poly data
POINT_DENSITY = 100;
r_lim = [min(r),max(r)];
r_poly = linspace(r_lim(1),r_lim(2),POINT_DENSITY);
x_poly = Disp_Poly.evaluate_polynomial(r_poly);
x_poly = x_poly(outputs,:);
f_poly = Force_Poly.evaluate_polynomial(r_poly);
h0_poly = zeros(1,POINT_DENSITY);

%h  poly data
H_SCALE_FACTOR = 1;
h_poly = r_poly*H_SCALE_FACTOR;
[R,H] = meshgrid(r_poly,h_poly);
num_points = POINT_DENSITY^2;
H0 = zeros(POINT_DENSITY);


R_linear = reshape(R,1,num_points);
H_linear = reshape(H,1,num_points);
H0_linear = reshape(H0,1,num_points);

H_DISP_linear = H_Displacement_Poly.evaluate_polynomial([R_linear;H0_linear;H_linear]);



figure
tiledlayout("flow")
for iTile = 1:NUM_OUTPUTS
    
    H_DISP = reshape(H_DISP_linear(outputs(iTile),:),size(H));
    nexttile
    box on
    hold on
    plot3(r,h,x(iTile,:),'kx')
    plot3(r_poly,h0_poly,x_poly(iTile,:))
    plot3(r,h_L_displacement,perturbation(outputs(iTile),:) + x(iTile,:),'rx')
    mesh(R,H,H_DISP)
    hold off
    xlabel("r")
    ylabel("h")
    zlabel("x_{" + outputs(iTile) + "}")
end


H_FORCE_linear = H_Force_Poly.evaluate_polynomial([R_linear;H0_linear;H_linear]);

figure
tiledlayout("flow")
for iTile = 1:3
    nexttile
    H_FORCE = reshape(H_FORCE_linear(iTile,:),size(H));
    if iTile > 1
        f = zeros(size(r));
        f_poly = zeros(size(r_poly));
    end
    box on
    hold on
    plot3(r,h,f,'kx')
    plot3(r_poly,h0_poly,f_poly)
    plot3(r,h_L_displacement,f_L,'rx')
    mesh(R,H,H_FORCE)
    hold off
    xlabel("r")
    ylabel("h")
    zlabel("f_" + iTile)
end
