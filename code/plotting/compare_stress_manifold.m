function ax = compare_stress_manifold(manifolds,varargin)
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);
Plot_Settings = struct([]);
ax = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "opts"
            Plot_Settings = keyword_values{arg_counter};
        case "axes"
            ax = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
if isempty(ax)
    fig = figure;
    ax = axes(fig);
end
if ~isfield(Plot_Settings,"light_on")
    Plot_Settings.light_on = 0;
end
if ~isfield(Plot_Settings,"energy_limit")
    Plot_Settings.energy_limit = 0;
end

if ~isfield(Plot_Settings,"legend")
    Plot_Settings.legend = 1;
end

if ~isfield(Plot_Settings,"mesh_alpha")
    Plot_Settings.mesh_alpha = 1;
end


colours = get_plot_colours(1);

num_manifolds = length(manifolds);
for iManifold = 1:num_manifolds
    Manifold = manifolds{iManifold};
    Dyn_Data = Manifold.system;
    validation_manifold_plotting = false;
    if isfield(Manifold,"plot_validation_manifold")
        validation_manifold_plotting = Manifold.plot_validation_manifold;
    end
    validation_orbit_plotting = true;
    if isfield(Manifold,"plot_validation_orbit")
        validation_orbit_plotting = Manifold.plot_validation_orbit;
    end

    if isstring(Dyn_Data)
        Dyn_Data = initalise_dynamic_data(Dyn_Data);
    end
    ax = plot_manifold(ax,Dyn_Data,Plot_Settings,colours(1,:));

    if ~isfield(Manifold,"orbit")
        continue
    end

    
    num_orbits = size(Manifold.orbit,1);
    for iOrbit = 1:num_orbits
        orbit_data = Manifold.orbit(iOrbit,:);
        if isstring(orbit_data(1))
            switch orbit_data(2)
                case "closest"
                    manifold_id = str2double(orbit_data(3));
                    sol_num = str2double(orbit_data(1));
                    comparative_manifold = manifolds{manifold_id};
                    Dyn_Data_Comp = comparative_manifold.system;
                    comparative_orbit = comparative_manifold.orbit(iOrbit,:);
                    orbit_data = get_closest_orbit(Dyn_Data,sol_num,Dyn_Data_Comp,comparative_orbit);
                otherwise
                    continue
            end
        end
        for jOrbit = 1:size(orbit_data,1)
            ax = plot_orbit(ax,Dyn_Data,orbit_data(jOrbit,:),Plot_Settings,validation_orbit_plotting,validation_manifold_plotting);
        end
    end

end


switch Plot_Settings.plot_type
    case "physical"
        label_base = "x_";
    case "modal"
        label_base = "q_";
    otherwise
        error("Unknown plot type: " + Plot_Settings.plot_type)
end

labels = arrayfun(@(iLabel) label_base + "{" + iLabel + "}",Plot_Settings.coords);

box(ax,"on")
xlabel(ax,labels(1));
ylabel(ax,labels(2))
zlabel(ax,labels(3))




if Plot_Settings.light_on
    light(ax,"Position",[0,0,1],"Color",0.8*[1,1,1])
end
% daspect([1 1 1])
%-------------------------------------------------------------------------%
if ~Plot_Settings.legend
    return %#ok<*UNRCH>
end
lines = ax.Children;
num_lines = size(lines,1);
for iLine = 1:num_lines
    line = lines(iLine);

    line_tag = line.Tag;
    if isempty(line_tag)
        if isprop(line,"Annotation")
            line.Annotation.LegendInformation.IconDisplayStyle = "off";
        end
        continue
    end
    tag_id = split(line_tag,"-");
    switch tag_id{1}
        case {"outline","grid_line"}
            line.HandleVisibility = "off";
            continue
    end
    manifold_def = "_{" + tag_id{2} + "}";
    switch tag_id{1}
        case "m"
            line.DisplayName = "$\mathcal{R}" + manifold_def + "$";
        case {"vo","o"}
            orbit_def = "[" + join(string(tag_id{3}),",") + "]";
            if tag_id{1} == "o"
                orbit_display = "\tilde{o}";
            elseif tag_id{1} == "vo"
                orbit_display = "\hat{o}";
            end
            line.DisplayName = "$" + orbit_display + manifold_def + ":" + orbit_def + "$";
        case "vm"
            line.DisplayName = "$\mathcal{V}(t)" + manifold_def + "$";
        
    end
end
leg = legend(ax);
leg.Interpreter = "latex";
end
%-------------------------------------------------------------------------%
function ax = plot_manifold(ax,Dyn_Data,Plot_Settings,colour)
PLOT_MODE = "ellipse";
mesh_alpha = Plot_Settings.mesh_alpha;
LINE_WIDTH = 2;
PLOT_RESOLUTION = 21;


mesh_settings = {"EdgeColor","none","EdgeAlpha",0,"LineWidth",1,"FaceColor","interp","FaceLighting","gouraud","FaceAlpha",mesh_alpha};
line_settings = {"k-","linewidth",LINE_WIDTH};

Rom = Dyn_Data.Dynamic_Model;

manifold_name = "m-" + join(string(Rom.Model.reduced_modes),",");

num_r_modes = length(Rom.Model.reduced_modes);
plot_order = Plot_Settings.coords;

hold(ax,"on")
switch num_r_modes
    case 1
        r_lim = Rom.reduced_displacement_limits;
        r = linspace(r_lim(1),r_lim(2),PLOT_RESOLUTION);
        if Plot_Settings.energy_limit ~= 0
            potential = Rom.Potential_Polynomial.evaluate_polynomial(r);
            energy_lim = Plot_Settings.energy_limit;
            positive_limit = interp1(potential(sign(r) == 1),r(sign(r) == 1),energy_lim);
            negative_limit = interp1(potential(sign(r) == -1),r(sign(r) == -1),energy_lim);
            r(potential > energy_lim) = [];
            r = [negative_limit,r,positive_limit];
        end
        x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r);


        plot3(ax,x_tilde(plot_order(1),:),x_tilde(plot_order(2),:),x_tilde(plot_order(3),:),line_settings{:},"tag",manifold_name);
    case 2
        switch PLOT_MODE
            case "rectangle"
                Z_BC = rectangle_manifold_mesh(Rom,Plot_Settings);
                z_grid = Z_BC(:,:,plot_order);
            case "ellipse"
                Phy_Poly = Rom.Physical_Displacement_Polynomial;
                Potential_Poly = Rom.Potential_Polynomial;
                energy_lim = Plot_Settings.energy_limit;
                [x_grid,y_grid] = get_2d_plotting_data(Phy_Poly,Potential_Poly,energy_lim);
                num_points = numel(x_grid);
                x_lin = reshape(x_grid,1,num_points);
                y_lin = reshape(y_grid,1,num_points);
                z_lin = Phy_Poly.evaluate_polynomial([x_lin;y_lin],plot_order);
                z_grid = reshape(z_lin',size(x_grid,1),size(x_grid,2),3);
        end
        colour_data = ones([size(z_grid,[1,2]),3]).*reshape(colour,[1,1,3]);
        mesh(ax,z_grid(:,:,1),z_grid(:,:,2),z_grid(:,:,3),colour_data,mesh_settings{:},"tag",manifold_name);

        switch PLOT_MODE
            case "ellipse"
                grid_line_style = {"Color",[0.25,0.25,0.25,0.2],"LineWidth",0.2,"Tag","grid_line"};
                outline_style = {"LineWidth",1,"Tag","outline","Color","k"};
                num_radii = size(z_grid,2);
                for iRadius = 1:(num_radii-1)
                    plot3(ax,z_grid(:,iRadius,1),z_grid(:,iRadius,2),z_grid(:,iRadius,3),grid_line_style{:})
                end
                plot3(ax,z_grid(:,end,1),z_grid(:,end,2),z_grid(:,end,3),outline_style{:})

        end
end
hold(ax,"off")
end
%--


%--
function Z_BC = rectangle_manifold_mesh(Rom,Plot_Settings)
PLOT_RESOLUTION = 21;
num_points = PLOT_RESOLUTION^2;

limits = Rom.Physical_Displacement_Polynomial.input_limit;
poly_bound = polyshape(limits');
[x_lim,y_lim] = boundingbox(poly_bound);


x = linspace(x_lim(1),x_lim(2),PLOT_RESOLUTION);
y = linspace(y_lim(1),y_lim(2),PLOT_RESOLUTION);
[X,Y] = meshgrid(x,y);

if Plot_Settings.energy_limit
    X_BC = nan(PLOT_RESOLUTION);
    Y_BC = nan(PLOT_RESOLUTION);
    for iCol = 1:PLOT_RESOLUTION
        x_vec = [X(:,iCol),Y(:,iCol)];
        valid_point = isinterior(poly_bound,x_vec);
        X_BC(valid_point,iCol) = X(valid_point,iCol);
        Y_BC(valid_point,iCol) = Y(valid_point,iCol);
    end
else
    X_BC = X;
    Y_BC = Y;
end

X_array = reshape(X_BC,1,num_points);
Y_array = reshape(Y_BC,1,num_points);
Z_array = Rom.Physical_Displacement_Polynomial.evaluate_polynomial([X_array;Y_array]);
num_dof = size(Z_array,1);
Z_BC = reshape(Z_array,[num_dof,size(X_BC)]);
Z_BC = permute(Z_BC,[2,3,1]);
end
%-------------------------------------------------------------------------%
function ax = plot_orbit(ax,Dyn_Data,orbit_data,Plot_Settings,validation_orbit_plotting,validation_manifold_plotting)
LINE_WIDTH = 2;

line_style = {"linewidth",LINE_WIDTH};
orbit_style = "-";
Rom = Dyn_Data.Dynamic_Model;


solution_num = orbit_data(1);
orbit_id = orbit_data(2);
[Orbit, Validation_orbit] = Dyn_Data.get_orbit(solution_num,orbit_id,1);

orbit_name = "o-" + join(string(Rom.Model.reduced_modes),",") + "-" + solution_num + "," + orbit_id;

num_r_modes = length(Rom.Model.reduced_modes);
if num_r_modes == 1
    orbit_style = ".-";
end
r_orbit = Orbit.xbp(:,1:num_r_modes)';
t0 = Orbit.tbp';

x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);

plot_order = Plot_Settings.coords;
hold(ax,"on")
p = plot3(ax,x_tilde(plot_order(1),:),x_tilde(plot_order(2),:),x_tilde(plot_order(3),:),orbit_style,line_style{:},"tag",orbit_name);
hold(ax,"off")
data_tip_row = dataTipTextRow("time (s)",t0);
p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

if isempty(Validation_orbit) || ~validation_orbit_plotting
    return
end
validation_orbit_name = "v" + orbit_name;
Static_Data = load_static_data(Rom);
Sol = Dyn_Data.load_solution(solution_num,"validation");
Static_Data = Static_Data.add_validation_data(Sol.validation_modes);
Rom = Reduced_System(Static_Data);

dof = size(x_tilde,1);
h_orbit = Validation_orbit.h;
num_points = size(h_orbit,2);
x_hat = zeros(dof,num_points);
Theta_Hat_Grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial;
for iPoint = 1:num_points
    x_hat_grad_orbit = Theta_Hat_Grad.evaluate_polynomial(r_orbit(:,iPoint));
    x_hat(:,iPoint) = x_tilde(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
end

hold(ax,"on")
p = plot3(ax,x_hat(plot_order(1),:),x_hat(plot_order(2),:),x_hat(plot_order(3),:),"-",line_style{:},"tag",validation_orbit_name);
hold(ax,"off")
p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

if validation_manifold_plotting > 0
    % Static_Data = load_static_data(Rom);
    % Static_Data = Static_Data.add_validation_data(Sol.validation_modes,1);
    % Rom = Reduced_System(Static_Data);
    t_point = validation_manifold_plotting;
    ax = plot_validation_manifold_time(ax,t_point,Rom,r_orbit,h_orbit,Plot_Settings);
end
end
%-------------------------------------------------------------------------%
function ax = plot_validation_manifold_time(ax,t_point,Rom,r_orbit,h_orbit,Plot_Settings)
H_LIM_SCALE_FACTOR = 1.5;
PLOT_RESOLUTION = 2;
MESH_ALPHA = 1;

manifold_colour = get_plot_colours(3);

mesh_settings = {"EdgeColor","none","FaceColor","interp","FaceLighting","gouraud","FaceAlpha",MESH_ALPHA};
manifold_name = "vm-" + join(string(Rom.Model.reduced_modes),",");

plot_order = Plot_Settings.coords;

num_h_modes = size(h_orbit,1);
if num_h_modes > 2
    return
end
%---------------------------------------
orbit_lim = [min(h_orbit,[],2),max(h_orbit,[],2)].*H_LIM_SCALE_FACTOR;

num_points = PLOT_RESOLUTION^2;
x = linspace(orbit_lim(1,1),orbit_lim(1,2),PLOT_RESOLUTION);
y = linspace(orbit_lim(2,1),orbit_lim(2,2),PLOT_RESOLUTION);
[X,Y] = meshgrid(x,y);

X_array = reshape(X,1,num_points);
Y_array = reshape(Y,1,num_points);
h_array = [X_array;Y_array];

x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit(:,t_point));
x_hat_grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_orbit(:,t_point));

x_hat = x_tilde + x_hat_grad*h_array;

num_dof = size(x_hat,1);
Z = reshape(x_hat,[num_dof,size(X)]);
Z = permute(Z,[2,3,1]);

colour_data = ones([size(Z,[1,2]),3]).*reshape(manifold_colour,[1,1,3]);

%validation whisker
hold(ax,"on")
mesh(ax,Z(:,:,plot_order(1)),Z(:,:,plot_order(2)),Z(:,:,plot_order(3)),colour_data,mesh_settings{:},"tag",manifold_name);
hold(ax,"off")

%r orbit time point
% hold(ax,"on")
% plot3(ax,x_tilde(plot_order(1)),x_tilde(plot_order(2)),x_tilde(plot_order(3)),"x")
% hold(ax,"off")
% 
% %h orbit time point
% h_t = h_orbit(:,t_point);
% x_hat_t = x_tilde + x_hat_grad*h_t;
% hold(ax,"on")
% plot3(x_hat_t(plot_order(1)),x_hat_t(plot_order(2)),x_hat_t(plot_order(3)),"x")
% hold(ax,"off")
end



function ax = plot_validation_manifold_orbits(ax,Rom,Dyn_Data,solution_num,Plot_Settings)
Sol = Dyn_Data.load_solution(solution_num);
num_orbits = Sol.num_orbits;
num_r_modes = length(Rom.Model.reduced_modes);
plot_order = Plot_Settings.coords;

hold(ax,"on")
for iOrbit = 1:num_orbits
    [Orbit, Validation_Orbit] = Dyn_Data.get_orbit(solution_num,iOrbit,1);

    r_orbit = Orbit.xbp(:,1:num_r_modes)';
    x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);
    dof = size(x_tilde,1);
    h_orbit = Validation_Orbit.h;
    num_points = size(h_orbit,2);
    x_hat = zeros(dof,num_points);
    Theta_Hat_Grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial;
    for iPoint = 1:num_points
        x_hat_grad_orbit = Theta_Hat_Grad.evaluate_polynomial(r_orbit(:,iPoint));
        x_hat(:,iPoint) = x_tilde(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
    end

    
    plot3(ax,x_hat(plot_order(1),:),x_hat(plot_order(2),:),x_hat(plot_order(3),:),"-","LineWidth",0.5,"Color",[0.3,0.3,0.3]);
   
end
hold(ax,"off")
end
% function ax = plot_validation_manifold(ax,Rom,orbit_data,r_orbit,h_orbit,Plot_Settings)

% NUM_DIMENSIONS = 2;
% H_LIM_SCALE_FACTOR = 1.5;
% PLOT_RESOLUTION = 101;
% MESH_ALPHA = 0.5;
%
% mesh_settings = {"EdgeColor","none","FaceColor","interp","FaceLighting","gouraud","FaceAlpha",MESH_ALPHA};
% manifold_name = "vm-" + join(string(Rom.Model.reduced_modes),",");
%
% plot_order = Plot_Settings.coords;
% %---------------------------------------
% r_transform = Rom.Model.reduced_eigenvectors'*Rom.Model.mass;
%
% switch NUM_DIMENSIONS
%     case 1
%         orbit_lim = [min(h_orbit(2,:)),max(h_orbit(2,:))].*H_LIM_SCALE_FACTOR;
%         x = linspace(orbit_lim(1),orbit_lim(2),PLOT_RESOLUTION);
%
%
%         num_t_points = size(r_orbit,2);
%         hold(ax,"on")
%         for iPoint = 1:num_t_points
%             r_i = r_orbit(:,iPoint);
%
%             theta_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_i);
%             displacement_gradient = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_i);
%             target_h = h_orbit(:,iPoint);
%             theta_hat = theta_tilde + displacement_gradient*target_h;
%
%             theta_plot = [theta_tilde,theta_hat];
%             theta_grad = diff(theta_plot,1,2);
%             theta_line = theta_tilde + theta_grad.*[-H_LIM_SCALE_FACTOR,H_LIM_SCALE_FACTOR];
%
%             plot3(ax,theta_line(plot_order(1),:),theta_line(plot_order(2),:),theta_line(plot_order(3),:));
%         end
%         hold(ax,"off")
%     case 2
%         orbit_lim = [min(h_orbit,[],2),max(h_orbit,[],2)].*H_LIM_SCALE_FACTOR;
%
%         num_points = PLOT_RESOLUTION^2;
%         x = linspace(orbit_lim(1,1),orbit_lim(1,2),PLOT_RESOLUTION);
%         y = linspace(orbit_lim(2,1),orbit_lim(2,2),PLOT_RESOLUTION);
%         [X,Y] = meshgrid(x,y);
%
%         X_array = reshape(X,1,num_points);
%         Y_array = reshape(Y,1,num_points);
%         h_array = [X_array;Y_array];
%
%
%         num_t_points = size(r_orbit,2);
%
%         Validation_Input = Rom.get_solver_inputs("h_prediction");
%
%         num_r_modes = size(r_orbit,1)/2;
%         r = r_orbit(1:num_r_modes,:);
%         r_dot = r_orbit((num_r_modes+1):end,:);
%
%         Eom_Input = Rom.get_solver_inputs("coco_backbone");
%
%
%
%         hold(ax,"on")
%         for iPoint = 1:num_t_points
%             r_i = r(:,iPoint);
%             r_dot_i = r_dot(:,iPoint);
%
%             z_ddot = coco_eom(0,r_orbit(:,iPoint),0,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
%             r_ddot_i = z_ddot((num_r_modes+1):end,:);
%             [h_inertia,h_conv,h_stiff,h_force] = get_h_error_terms(r_i,r_dot_i,r_ddot_i,Validation_Input);
%
%
%             h = h_stiff\h_force;
%
%             displacement_gradient = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_i);
%             theta_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_i);
%             Z_array = theta_tilde + displacement_gradient*h;
%             plot3(Z_array(plot_order(1),:),Z_array(plot_order(2),:),Z_array(plot_order(3),:),"x")
%             % num_dof = size(Z_array,1);
%             % Z = reshape(Z_array,[num_dof,size(X)]);
%             % Z = permute(Z,[2,3,1]);
%             % colour_data = ones(size(Z,[1,2]));
%             % mesh(ax,Z(:,:,plot_order(1)),Z(:,:,plot_order(2)),Z(:,:,plot_order(3)),colour_data,mesh_settings{:},"tag",manifold_name);
%         end
%         hold(ax,"off")
% end

% end


function matched_orbits = get_closest_orbit(Dyn_Data,sol_num,Dyn_Data_Comp,comparative_orbit)
if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
if isstring(Dyn_Data_Comp)
    Dyn_Data_Comp = initalise_dynamic_data(Dyn_Data_Comp);
end

comparative_sol_num = comparative_orbit(1);
comparative_orbit_id = comparative_orbit(2);
Sol_Comp = Dyn_Data_Comp.load_solution(comparative_sol_num);
comparative_orbit_frequency = Sol_Comp.frequency(comparative_orbit_id);

Sol = Dyn_Data.load_solution(sol_num);
frequency = Sol.frequency;

% find sign changes
frequency_diff = frequency - comparative_orbit_frequency;
frequency_diff_sign = sign(frequency_diff);
sign_diff = find(diff(frequency_diff_sign) ~= 0);
num_matches = length(sign_diff);
matched_orbits = zeros(num_matches,2);
matched_orbits(:,1) = sol_num;
for iMatch = 1:num_matches
    matched_index = sign_diff(iMatch) + [0,1];
    bounds = frequency_diff(matched_index);
    [~,closest_bound] = min(abs(bounds));
    matched_orbits(iMatch,2) = matched_index(closest_bound);
end
end