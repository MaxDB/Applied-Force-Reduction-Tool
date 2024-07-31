clear
close all
Static_Data_1 = load("data\beam_oscillator_1\Static_Data.mat");
Static_Data_1 = Static_Data_1.Static_Data;

Static_Data_12 = load("data\beam_oscillator_12\Static_Data.mat");
Static_Data_12 = Static_Data_12.Static_Data;
%-------------------------------------------------------------------------%
Static_Data_1 = Static_Data_1.add_validation_data(2);


Rom_1 = Reduced_System(Static_Data_1);
Rom_12 = Reduced_System(Static_Data_12);
Disp_Poly_1 = Rom_1.Physical_Displacement_Polynomial;
H_Displacement_Poly_1 = Rom_1.Low_Frequency_Displacement_Polynomial;
Force_Poly_1 = Rom_1.Force_Polynomial;
H_Force_Poly_1 = Rom_1.Low_Frequency_Force_Polynomial;
Disp_Poly_12 = Rom_12.Physical_Displacement_Polynomial;

Model = Rom_12.Model;
%static data
r_1 = Static_Data_1.reduced_displacement;
% h_1 = zeros(size(r_1));
% x_1 = Static_Data_1.physical_displacement;
% f_1 = Static_Data_1.restoring_force;

%poly data
POINT_DENSITY = 100;
limits = Rom_12.Physical_Displacement_Polynomial.input_limit;
poly_bound = polyshape(limits');
[x_lim,y_lim] = boundingbox(poly_bound);


x = linspace(x_lim(1),x_lim(2),POINT_DENSITY);
y = linspace(y_lim(1),y_lim(2),POINT_DENSITY);
[R1_12,R2_12] = meshgrid(x,y);

R1_12_BC = nan(POINT_DENSITY);
R2_12_BC = nan(POINT_DENSITY);
for iCol = 1:POINT_DENSITY
    x_vec = [R1_12(:,iCol),R2_12(:,iCol)];
    valid_point = isinterior(poly_bound,x_vec);
    R1_12_BC(valid_point,iCol) = R1_12(valid_point,iCol);
    R2_12_BC(valid_point,iCol) = R2_12(valid_point,iCol);
end

r_1_lim = [min(r_1),max(r_1)];
r_poly_1 = linspace(r_1_lim(1),r_1_lim(2),POINT_DENSITY);

x_poly_1 = Disp_Poly_1.evaluate_polynomial(r_poly_1);
s_poly_1 = Model.reduced_eigenvectors(:,2)'*Model.mass*x_poly_1;
f_poly_1 = Force_Poly_1.evaluate_polynomial(r_poly_1);
h0_poly_1 = zeros(1,POINT_DENSITY);

%h  poly data
num_points = POINT_DENSITY^2;
H0 = zeros(POINT_DENSITY);


R_1_linear = reshape(R1_12_BC,1,num_points);
S_1_linear = Model.reduced_eigenvectors(:,2)'*Model.mass*Disp_Poly_1.evaluate_polynomial(R_1_linear);
R_2_linear = reshape(R2_12_BC,1,num_points);
H0_linear_1 = reshape(H0,1,num_points);

H_DISP_linear_1 = H_Displacement_Poly_1.evaluate_polynomial([R_1_linear;H0_linear_1;R_2_linear]);
DISP_linear_12 = Disp_Poly_12.evaluate_polynomial([R_1_linear;R_2_linear]);

S_1 = reshape(S_1_linear,size(R1_12_BC));

figure
tiledlayout("flow")
for iDof = 1:2
    nexttile
    H_DISP_1 = reshape(H_DISP_linear_1(iDof,:),size(R1_12_BC));
    DISP_12 = reshape(DISP_linear_12(iDof,:),size(R1_12_BC));
    box on
    hold on
  
    plot3(r_poly_1,s_poly_1,x_poly_1(iDof,:))
    
    surf(R1_12_BC,R2_12_BC+S_1,H_DISP_1,ones(POINT_DENSITY),"EdgeColor",get_plot_colours("grey"))
    surf(R1_12_BC,R2_12_BC,DISP_12,ones(POINT_DENSITY) + 1,"EdgeColor",get_plot_colours("grey"))
    colormap(get_plot_colours(1:2))
    hold off
    xlabel("q_1")
    ylabel("q_2")
    zlabel("x_"+iDof)
end

Potential_Poly_1 = Rom_1.Potential_Polynomial;
H_Potential_Poly_1 = integrate_polynomial(Rom_1.Low_Frequency_Force_Polynomial);
Potential_Poly_12 = Rom_12.Potential_Polynomial;

v_poly_1 = Potential_Poly_1.evaluate_polynomial(r_poly_1);

R1_12_BC = nan(POINT_DENSITY);
R2_12_BC = nan(POINT_DENSITY);
for iCol = 1:POINT_DENSITY
    x_vec = [R1_12(:,iCol),R2_12(:,iCol)];
    valid_point = Potential_Poly_12.evaluate_polynomial(x_vec') < Rom_12.Model.energy_limit;
    R1_12_BC(valid_point,iCol) = R1_12(valid_point,iCol);
    R2_12_BC(valid_point,iCol) = R2_12(valid_point,iCol);
end

R_1_linear = reshape(R1_12_BC,1,num_points);
R_2_linear = reshape(R2_12_BC,1,num_points);
H0_linear_1 = reshape(H0,1,num_points);
S_1_linear = Model.reduced_eigenvectors(:,2)'*Model.mass*Disp_Poly_1.evaluate_polynomial(R_1_linear);

H_POTENTIAL_linear_1 = H_Potential_Poly_1.evaluate_polynomial([R_1_linear;H0_linear_1;R_2_linear]);
POTENTIAL_linear_12 = Potential_Poly_12.evaluate_polynomial([R_1_linear;R_2_linear]);

H_POTENTIAL_1 = reshape(H_POTENTIAL_linear_1,size(R1_12_BC));
POTENTIAL_12 = reshape(POTENTIAL_linear_12,size(R1_12_BC));
S_1 = reshape(S_1_linear,size(R1_12_BC));

figure
box on 
hold on
plot3(r_poly_1,h0_poly_1+s_poly_1,v_poly_1)
surf(R1_12_BC,R2_12_BC+S_1,H_POTENTIAL_1,ones(POINT_DENSITY),'EdgeColor',get_plot_colours("grey"))
surf(R1_12_BC,R2_12_BC,POTENTIAL_12,ones(POINT_DENSITY)+1,'EdgeColor',get_plot_colours("grey"))
colormap(get_plot_colours(1:2))
hold off
xlabel("q_1")
ylabel("q_2")
zlabel("V")

% H_FORCE_linear = H_Force_Poly_1.evaluate_polynomial([R_linear_1;H0_linear_1;H_linear_1]);
% H_FORCE = reshape(H_FORCE_linear(1,:),size(H_1));
% 
% figure
% box on
% hold on
% plot3(r_poly_1,h0_poly_1,f_poly_1)
% mesh(R_1,H_1,H_FORCE)
% hold off
