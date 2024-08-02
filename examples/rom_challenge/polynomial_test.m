clear
close all
load("data\exhaust_17\Static_Data\Static_Data.mat");
%-------------------------------------------------------------------------%
Static_Data = Static_Data.add_validation_data(2);


plot_static_data("displacement",Static_Data)


Rom = Reduced_System(Static_Data);
Disp_Poly = Rom.Physical_Displacement_Polynomial;
H_Displacement_Poly = Rom.Low_Frequency_Displacement_Polynomial;
Force_Poly = Rom.Force_Polynomial;
H_Force_Poly = Rom.Low_Frequency_Force_Polynomial;


%static data
X_DOF = 2;
r = Static_Data.reduced_displacement;
h = zeros(size(r));
x = Static_Data.physical_displacement(X_DOF,:);
f = Static_Data.restoring_force;

%poly data
POINT_DENSITY = 100;
r_lim = [min(r),max(r)];
r_poly = linspace(r_lim(1),r_lim(2),POINT_DENSITY);
x_poly = Disp_Poly.evaluate_polynomial(r_poly);
x_poly = x_poly(X_DOF,:);
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
H_DISP = reshape(H_DISP_linear(X_DOF,:),size(H));


figure
box on
hold on
plot3(r,h,x,'x')
plot3(r_poly,h0_poly,x_poly)
mesh(R,H,H_DISP)
hold off


H_FORCE_linear = H_Force_Poly.evaluate_polynomial([R_linear;H0_linear;H_linear]);
H_FORCE = reshape(H_FORCE_linear(1,:),size(H));

figure
box on
hold on
plot3(r,h,f,'x')
plot3(r_poly,h0_poly,f_poly)
mesh(R,H,H_FORCE)
hold off
