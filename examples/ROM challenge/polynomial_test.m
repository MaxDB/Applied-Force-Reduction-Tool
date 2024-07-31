clear
close all
load("data\exhaust_1\Static_Data.mat");
%-------------------------------------------------------------------------%
H_MODE = 7;

Static_Data = Static_Data.add_validation_data(H_MODE);
Rom = Reduced_System(Static_Data);

Model = Rom.Model;

h_mode_transform = [Model.reduced_eigenvectors,Model.low_frequency_eigenvectors(:,Model.low_frequency_modes == H_MODE)]'*Model.mass;

Disp_Poly = Rom.Physical_Displacement_Polynomial;
H_Displacement_Poly = Rom.Low_Frequency_Displacement_Polynomial;
Force_Poly = Rom.Force_Polynomial;
H_Force_Poly = Rom.Low_Frequency_Force_Polynomial;


%static data
NUM_OUTPUTS = 4;
outputs = randi(size(Disp_Poly,1),1,NUM_OUTPUTS);

r = Static_Data.reduced_displacement;
h = zeros(size(r));
x = Static_Data.physical_displacement;
s = h_mode_transform(2,:)*x;
x = x(outputs,:);
f = Static_Data.restoring_force;

perturbation_1 = squeeze(Static_Data.condensed_perturbation(:,1,:));
perturbation_2 = squeeze(Static_Data.condensed_perturbation(:,2,:));

h_displacement = Static_Data.get_h_displacement;
h_r_1 = squeeze(h_displacement(1,1,:))';
h_L_1 = squeeze(h_displacement(2,1,:))';
h_r_2 = squeeze(h_displacement(1,2,:))';
h_L_2 = squeeze(h_displacement(2,2,:))';


f_L = Static_Data.perturbation_scale_factor*ones(1,size(f,2));


%poly data
POINT_DENSITY = 100;
r_lim = [min(r),max(r)];
r_poly = linspace(r_lim(1),r_lim(2),POINT_DENSITY);
x_poly = Disp_Poly.evaluate_polynomial(r_poly);
s_poly = h_mode_transform(2,:)*x_poly;
x_poly = x_poly(outputs,:);
f_poly = Force_Poly.evaluate_polynomial(r_poly);
h0_poly = zeros(1,POINT_DENSITY);

%h  poly data
H_SCALE_FACTOR = 0.1;
h_poly = r_poly*H_SCALE_FACTOR;
[R,H] = meshgrid(r_poly,h_poly);
num_points = POINT_DENSITY^2;
H0 = zeros(POINT_DENSITY);


R_linear = reshape(R,1,num_points);
H_linear = reshape(H,1,num_points);
H0_linear = reshape(H0,1,num_points);

H_DISP_linear = H_Displacement_Poly.evaluate_polynomial([R_linear;H0_linear;H_linear]);
S_linear = h_mode_transform(2,:)*Disp_Poly.evaluate_polynomial(R_linear);
S = reshape(S_linear,size(H));

figure
tiledlayout("flow")
for iTile = 1:NUM_OUTPUTS
    
    H_DISP = reshape(H_DISP_linear(outputs(iTile),:),size(H));
    nexttile
    box on
    hold on
    plot3(r,h+s,x(iTile,:),'kx')
    plot3(r_poly,h0_poly + s_poly,x_poly(iTile,:))
    % plot3(r+h_r_1,h_L_1 + s,perturbation_1(outputs(iTile),:) + x(iTile,:),'rx')
    plot3(r+h_r_2,h_L_2 + s,perturbation_2(outputs(iTile),:) + x(iTile,:),'rx')
    mesh(R,H+S,H_DISP)
    hold off
    xlabel("q_1")
    ylabel("q_{" + H_MODE + "}")
    zlabel("x_{" + outputs(iTile) + "}")
end


H_FORCE_linear = H_Force_Poly.evaluate_polynomial([R_linear;H0_linear;H_linear]);

f_r = [f;zeros(2,size(f,2))];
f_r_poly = [f_poly;zeros(2,POINT_DENSITY)];
f_h = [f;zeros(1,size(f,2));f_L];
figure
tiledlayout("flow")
for iTile = 1:3
    nexttile
    H_FORCE = reshape(H_FORCE_linear(iTile,:),size(H));

    box on
    hold on
    plot3(r,h+s,f_r(iTile,:),'kx')
    plot3(r_poly,h0_poly+s_poly,f_r_poly(iTile,:))
    % plot3(r + h_r_1,h_L_1 + s,f_L,'rx')
    plot3(r + h_r_2,h_L_2 + s,f_h(iTile,:),'rx')
    mesh(R,H+S,H_FORCE)
    hold off
    xlabel("q_1")
    ylabel("q_{" + H_MODE + "}")
    zlabel("f_" + iTile)
end
