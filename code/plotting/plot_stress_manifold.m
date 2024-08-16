function plot_stress_manifold(Dyn_Data,L_modes,varargin)
PLOT_H_R = 0;
PLOT_RESOLUTION = 51;
LIM_SF = 1;
PHYSICAL_COORDINATES = 1;
PLOT_MESH = 0;

%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

Dyn_Data_Comp = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "comparison"
            Dyn_Data_Comp = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%

if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
if isstring(Dyn_Data_Comp)
    Dyn_Data_Comp = initalise_dynamic_data(Dyn_Data_Comp);
end

Rom = Dyn_Data.Dynamic_Model;

Static_Data = load_static_data(Rom);
Static_Data = Static_Data.add_validation_data(L_modes);


Rom = Reduced_System(Static_Data);


if ~PHYSICAL_COORDINATES
    mass = Rom.Model.mass;
    stiffness = Rom.Model.stiffness;
    [evec,eval] = eig(stiffness,mass);
end

%get all stress manifold parameters
%generate cloud of physical displacement points

r_lim = Rom.reduced_displacement_limits;
h_r_lim = r_lim*LIM_SF;

if isempty(Dyn_Data_Comp)
    h_L_lim = r_lim*LIM_SF;
else
    h_L_lim = Dyn_Data_Comp.Dynamic_Model.reduced_displacement_limits(2,:)*LIM_SF;
end

if ~PLOT_H_R
    h_r_resolution = 1;
    h_r_lim = [0,0];
else
    h_r_resolution = PLOT_RESOLUTION; %#ok<*UNRCH>
end

r = linspace(r_lim(1),r_lim(2),PLOT_RESOLUTION);
h_r = linspace(h_r_lim(1),h_r_lim(2),PLOT_RESOLUTION);
h_L = linspace(h_L_lim(1),h_L_lim(2),PLOT_RESOLUTION);


x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r);
num_dofs = size(x_tilde,1);

figure
box on
hold on
if PHYSICAL_COORDINATES
    xlabel("x_1")
    ylabel("x_2")
    zlabel("x_3")
else
    xlabel("q_1")
    ylabel("q_2")
    zlabel("q_3")
end

if PLOT_MESH
    [R,H_L] = meshgrid(r,h_L);
    R_lin = reshape(R,PLOT_RESOLUTION^2,1);
    H_L_lin = reshape(H_L,PLOT_RESOLUTION^2,1);
   
    
    num_points = size(R_lin,1);
    X_HAT_lin = zeros(PLOT_RESOLUTION^2,num_dofs);
    for iPoint = 1:num_points
         x_hat_grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(R_lin(iPoint));
         X_HAT_lin(iPoint,:) = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(R_lin(iPoint)) + x_hat_grad*[0;H_L_lin(iPoint)];
    end

    X_HAT = reshape(X_HAT_lin,PLOT_RESOLUTION,PLOT_RESOLUTION,3);
    mesh(X_HAT(:,:,1),X_HAT(:,:,2),X_HAT(:,:,3))
else

    for iR = 1:PLOT_RESOLUTION
        x_hat_grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r(:,iR));
        for iH_R = 1:h_r_resolution
            x_hat = zeros(num_dofs,PLOT_RESOLUTION);
            for iH_L = 1:PLOT_RESOLUTION
                h = [h_r(:,iH_R);h_L(:,iH_L)];
                x_hat(:,iH_L) = x_tilde(:,iR) + x_hat_grad*h;
            end
            if ~PHYSICAL_COORDINATES
                x_hat = evec*x_hat;
            end

            switch num_dofs
                case 2
                    plot(x_hat(1,:),x_hat(2,:),'-')
                case 3
                    plot3(x_hat(1,:),x_hat(2,:),x_hat(3,:),'-')
            end

        end
    end
end
% x_hat_mesh = pc2surfacemesh(pointCloud(x_hat'),"poisson");
% surfaceMeshShow(x_hat_mesh)
if ~PHYSICAL_COORDINATES
    x_tilde = evec*x_tilde;
end
switch num_dofs
    case 2
        plot(x_tilde(1,:),x_tilde(2,:),'k-')
    case 3
        plot3(x_tilde(1,:),x_tilde(2,:),x_tilde(3,:),'k-')
end

hold off

PLOT_ORBIT = 0;
if PLOT_ORBIT
    orbit_id = 151;
    [orbit,validation_orbit] = Dyn_Data.get_orbit(2,orbit_id,1);

    r_orbit = orbit.xbp(:,1)';
    x_tilde_orbit = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);

    h_orbit = validation_orbit.h;
    num_points = size(h_orbit,2);
    x_hat_orbit = zeros(3,num_points);
    for iPoint = 1:num_points
        x_hat_grad_orbit = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_orbit(:,iPoint));
        x_hat_orbit(:,iPoint) = x_tilde_orbit(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
    end

    hold on
    plot3(x_tilde_orbit(1,:),x_tilde_orbit(2,:),x_tilde_orbit(3,:),'k-','LineWidth',2)
    plot3(x_hat_orbit(1,:),x_hat_orbit(2,:),x_hat_orbit(3,:),'r-','LineWidth',2)
hold off
end

if isempty(Dyn_Data_Comp)
    return
end

rom_2 = Dyn_Data_Comp.Dynamic_Model;

r_1 = r;
r_2_lim = rom_2.reduced_displacement_limits(2,:);
r_2 = linspace(r_2_lim(1),r_2_lim(2),PLOT_RESOLUTION);

hold on
for iR_1 = 1:PLOT_RESOLUTION
    r_1_plot = r_1(:,iR_1);
    x_tilde_12 = zeros(num_dofs,PLOT_RESOLUTION);
    for iR_2 = 1:PLOT_RESOLUTION
        r_2_plot = r_2(:,iR_2);
        x_tilde_12(:,iR_2) = rom_2.Physical_Displacement_Polynomial.evaluate_polynomial([r_1_plot;r_2_plot]);
    end
    if ~PHYSICAL_COORDINATES
        x_tilde_12 = evec*x_tilde_12;
    end
    switch num_dofs
        case 2
            plot(x_tilde_12(1,:),x_tilde_12(2,:),'-')
        case 3
            plot3(x_tilde_12(1,:),x_tilde_12(2,:),x_tilde_12(3,:),'-')
    end
end
hold off


end

