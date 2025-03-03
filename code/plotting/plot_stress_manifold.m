function ax = plot_stress_manifold(Dyn_Data,L_modes,varargin)
ALT_MODE = 0;
PLOT_H_R = 0;
PLOT_RESOLUTION = 51;
LIM_SF = 1;
PHYSICAL_COORDINATES = 1;

PLOT_LEGEND = 1;


LINE_WIDTH = 4;

MESH_ALPHA = 1;
MESH_COLOUR = 3;
COMPARISON_MESH_ALPHA = 0.6;
COMPARISON_MESH_COLOUR = 4;

%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

Dyn_Data_Comp = [];
orbit_id = [];
solution_num = [];
comparison_orbit_id = [];
comparison_solution_num = [];
ax = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "comparison"
            Dyn_Data_Comp = keyword_values{arg_counter};
        case "orbit"
            orbit_data = keyword_values{arg_counter};
            solution_num = orbit_data(1);
            orbit_id = orbit_data(2);
        case "comparison orbit"
            comparison_orbit_data = keyword_values{arg_counter};
            comparison_solution_num = comparison_orbit_data(1);
            comparison_orbit_id = comparison_orbit_data(2);
        case "axes"
            ax = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end

plot_mesh = isempty(Dyn_Data_Comp);
plot_mesh = 1;
plot_orbit = ~isempty(orbit_id);
plot_comparison_orbit = ~isempty(comparison_orbit_id);
if plot_orbit
    LIM_SF = 1;
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

if isempty(ax)
    figure
    box on
end
hold on
if PHYSICAL_COORDINATES
    xlabel("x_3")
    ylabel("x_2")
    zlabel("x_1")
else
    xlabel("q_1")
    ylabel("q_2")
    zlabel("q_3")
end



mesh_settings = {"EdgeColor","none","FaceColor","interp","FaceLighting","gouraud"};
if ~ALT_MODE
if plot_mesh
    [R,H_L] = meshgrid(r,h_L);
    R_lin = reshape(R,PLOT_RESOLUTION^2,1);
    H_L_lin = reshape(H_L,PLOT_RESOLUTION^2,1);
   
    
    num_points = size(R_lin,1);
    X_HAT_lin = zeros(PLOT_RESOLUTION^2,num_dofs);
    for iPoint = 1:num_points
         x_hat_grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(R_lin(iPoint));
         x_hat_lin = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(R_lin(iPoint)) + x_hat_grad*[0;H_L_lin(iPoint)];
         if ~PHYSICAL_COORDINATES
             x_hat_lin = evec'*mass*x_hat_lin;
         end
         X_HAT_lin(iPoint,:) = x_hat_lin;
    end

    X_HAT = reshape(X_HAT_lin,PLOT_RESOLUTION,PLOT_RESOLUTION,3);
    colour_data = zeros(size(X_HAT,[1,2]));
    mesh(X_HAT(:,:,3),X_HAT(:,:,2),X_HAT(:,:,1),colour_data,"FaceAlpha",MESH_ALPHA,mesh_settings{:},"Tag","validation")
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
                x_hat = evec'*mass*x_hat;
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
    x_tilde = evec'*mass*x_tilde;
end
switch num_dofs
    case 2
        plot(x_tilde(1,:),x_tilde(2,:),'k-',"LineWidth",LINE_WIDTH,"Tag","one mode")
    case 3
        plot3(x_tilde(3,:),x_tilde(2,:),x_tilde(1,:),'k-',"LineWidth",LINE_WIDTH,"Tag","one mode")
end

hold off


if plot_orbit
    [orbit,validation_orbit] = Dyn_Data.get_orbit(solution_num,orbit_id,1);

    r_orbit = orbit.xbp(:,1)';
    x_tilde_orbit = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);

    h_orbit = validation_orbit.h;
    num_points = size(h_orbit,2);
    x_hat_orbit = zeros(3,num_points);
    for iPoint = 1:num_points
        x_hat_grad_orbit = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_orbit(:,iPoint));
        x_hat_orbit(:,iPoint) = x_tilde_orbit(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
    end
    
    if ~PHYSICAL_COORDINATES
        x_tilde_orbit =evec'*mass*x_tilde_orbit;
        x_hat_orbit = evec'*mass*x_hat_orbit;
    end


    hold on
    plot3(x_tilde_orbit(3,:),x_tilde_orbit(2,:),x_tilde_orbit(1,:),'g-','LineWidth',2)
    plot3(x_hat_orbit(3,:),x_hat_orbit(2,:),x_hat_orbit(1,:),'r-','LineWidth',2)
    hold off
end



else
    BB_Sol = Dyn_Data.load_solution(solution_num);
    num_orbits = BB_Sol.num_orbits;
    num_r_modes = size(Static_Data,1);


    for iOrbit = 1:num_orbits
        [orbit,validation_orbit] = Dyn_Data.get_orbit(solution_num,iOrbit,1);
        r_orbit = orbit.xbp(:,1:num_r_modes)';
        x_tilde_orbit = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);

        h_orbit = validation_orbit.h;
        num_points = size(h_orbit,2);
        x_hat_orbit = zeros(3,num_points);
        for iPoint = 1:num_points
            x_hat_grad_orbit = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_orbit(:,iPoint));
            x_hat_orbit(:,iPoint) = x_tilde_orbit(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
        end

        if ~PHYSICAL_COORDINATES
            x_tilde_orbit =evec'*mass*x_tilde_orbit;
            x_hat_orbit = evec'*mass*x_hat_orbit;
        end


        plot3(x_tilde_orbit(1,:),x_tilde_orbit(2,:),x_tilde_orbit(3,:),'k.','LineWidth',2)
        plot3(x_hat_orbit(1,:),x_hat_orbit(2,:),x_hat_orbit(3,:),'r.')
    end
    hold off
end


light("Position",[1,1,1])

if isempty(Dyn_Data_Comp)
    return
end
PLOT_MESH = 1;


Rom_2 = Dyn_Data_Comp.Dynamic_Model;

r_1 = r;
r_2_lim = Rom_2.reduced_displacement_limits(2,:);
r_2 = linspace(r_2_lim(1),r_2_lim(2),PLOT_RESOLUTION);

if PLOT_MESH
    [R_1,R_2] = meshgrid(r_1,r_2);
    X_TILDE_12 = zeros(PLOT_RESOLUTION,PLOT_RESOLUTION,3);
    for iR_1 = 1:PLOT_RESOLUTION
        for iR_2 = 1:PLOT_RESOLUTION
            r_1_lin = R_1(iR_1,iR_2);
            r_2_lin = R_2(iR_1,iR_2);
            x_tilde_12_lin = Rom_2.Physical_Displacement_Polynomial.evaluate_polynomial([r_1_lin;r_2_lin]);
            
            if ~PHYSICAL_COORDINATES
                x_tilde_12_lin =evec'*mass*x_tilde_12_lin;
            end

            X_TILDE_12(iR_1,iR_2,:) = x_tilde_12_lin;
        end
    end
    
    
    colour_data = ones(size(X_TILDE_12,[1,2]));
    hold on
    mesh(X_TILDE_12(:,:,3),X_TILDE_12(:,:,2),X_TILDE_12(:,:,1),colour_data,"FaceAlpha",COMPARISON_MESH_ALPHA,mesh_settings{:},"Tag","two mode");
    hold off

else
    hold on
    for iR_1 = 1:PLOT_RESOLUTION
        r_1_plot = r_1(:,iR_1);
        x_tilde_12 = zeros(num_dofs,PLOT_RESOLUTION);
        for iR_2 = 1:PLOT_RESOLUTION
            r_2_plot = r_2(:,iR_2);
            x_tilde_12(:,iR_2) = Rom_2.Physical_Displacement_Polynomial.evaluate_polynomial([r_1_plot;r_2_plot]);
        end
        if ~PHYSICAL_COORDINATES
            x_tilde_12 = evec'*mass*x_tilde_12;
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


if plot_comparison_orbit
    orbit = Dyn_Data_Comp.get_orbit(comparison_solution_num,comparison_orbit_id);

    r_orbit = orbit.xbp(:,1)';
    x_tilde_orbit = Rom_2.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);

    % h_orbit = validation_orbit.h;
    % num_points = size(h_orbit,2);
    % x_hat_orbit = zeros(3,num_points);
    % for iPoint = 1:num_points
    %     x_hat_grad_orbit = Rom_2.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_orbit(:,iPoint));
    %     x_hat_orbit(:,iPoint) = x_tilde_orbit(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
    % end
    
    if ~PHYSICAL_COORDINATES
        x_tilde_orbit =evec'*mass*x_tilde_orbit;
        % x_hat_orbit = evec'*mass*x_hat_orbit;
    end


    hold on
    plot3(x_tilde_orbit(3,:),x_tilde_orbit(2,:),x_tilde_orbit(1,:),'b-','LineWidth',2)
    % plot3(x_hat_orbit(1,:),x_hat_orbit(2,:),x_hat_orbit(3,:),'r-','LineWidth',2)
    hold off
end


if ~isequal(MESH_COLOUR,0)
    colormap(get_plot_colours([MESH_COLOUR,COMPARISON_MESH_COLOUR]))
end

if PLOT_LEGEND
    ax = gca;
    lines = ax.Children;
    num_lines = size(lines,1);
    for iLine = 1:num_lines
        line = lines(iLine);
        if isa(line,"matlab.graphics.primitive.Light")
            continue
        end

        if line.Marker == "."
            line.DisplayName = "";
            continue
        end
        switch line.Tag
            case "one mode"
                % line.DisplayName = "$\mathcal W_{\{1\}}$";
                 line.DisplayName = "One mode stress manifold";
                uistack(line,"top")
                uistack(line,"down",2)
            case "two mode"
                % line.DisplayName = "$\mathcal W_{\{1,2\}}$";
                 line.DisplayName = "Two mode stress manifold";
                uistack(line,"top")
                
            case "validation"
                % line.DisplayName = "$\mathcal V_{\{1\}:\{1,2\}}$";
                line.DisplayName = "Two mode validation manifold";
                uistack(line,"top")
                uistack(line,"down",1)
        end
    end
    leg = legend;
    leg.Interpreter = "latex";
end
end
