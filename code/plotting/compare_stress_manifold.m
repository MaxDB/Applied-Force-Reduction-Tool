function ax = compare_stress_manifold(manifolds,varargin)
PLOT_LEGEND = 1;
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);
Plot_Settings = struct([]);

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "opts"
            Plot_Settings = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
fig = figure;
ax = axes(fig);


num_manifolds = length(manifolds);
for iManifold = 1:num_manifolds
    Manifold = manifolds{iManifold};
    Dyn_Data = Manifold.system;
    validation_manifold_plotting = false;
    if isfield(Manifold,"plot_validation_manifold")
        validation_manifold_plotting = Manifold.plot_validation_manifold;
    end

    if isstring(Dyn_Data)
        Dyn_Data = initalise_dynamic_data(Dyn_Data);
    end
    ax = plot_manifold(ax,Dyn_Data,Plot_Settings);

    if ~isfield(Manifold,"orbit")
        continue
    end
    num_orbits = size(Manifold.orbit,1);
    for iOrbit = 1:num_orbits
        orbit_data = Manifold.orbit(iOrbit,:);
        ax = plot_orbit(ax,Dyn_Data,orbit_data,Plot_Settings,validation_manifold_plotting);
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

light("Position",[1,1,1])
% daspect([1 1 1])
%-------------------------------------------------------------------------%
if ~PLOT_LEGEND
    return %#ok<*UNRCH>
end
lines = ax.Children;
num_lines = size(lines,1);
for iLine = 1:num_lines
    line = lines(iLine);

    line_tag = line.Tag;
    if isempty(line_tag)
        continue
    end
    tag_id = split(line_tag,"-");
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
    end
end
leg = legend;
leg.Interpreter = "latex";

%-------------------------------------------------------------------------%
    function ax = plot_manifold(ax,Dyn_Data,Plot_Settings)
        PLOT_RESOLUTION = 101;
        MESH_ALPHA = 0.5;
        LINE_WIDTH = 1;

        

        mesh_settings = {"EdgeColor","none","FaceColor","interp","FaceLighting","gouraud","FaceAlpha",MESH_ALPHA};
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

                x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r);
                
                
                plot3(ax,x_tilde(plot_order(1),:),x_tilde(plot_order(2),:),x_tilde(plot_order(3),:),line_settings{:},"tag",manifold_name);
            case 2
                num_points = PLOT_RESOLUTION^2;
                
                limits = Rom.Physical_Displacement_Polynomial.input_limit;
                poly_bound = polyshape(limits');
                [x_lim,y_lim] = boundingbox(poly_bound);


                x = linspace(x_lim(1),x_lim(2),PLOT_RESOLUTION);
                y = linspace(y_lim(1),y_lim(2),PLOT_RESOLUTION);
                [X,Y] = meshgrid(x,y);

                X_BC = nan(PLOT_RESOLUTION);
                Y_BC = nan(PLOT_RESOLUTION);
                for iCol = 1:PLOT_RESOLUTION
                    x_vec = [X(:,iCol),Y(:,iCol)];
                    valid_point = isinterior(poly_bound,x_vec);
                    X_BC(valid_point,iCol) = X(valid_point,iCol);
                    Y_BC(valid_point,iCol) = Y(valid_point,iCol);
                end

                X_array = reshape(X_BC,1,num_points);
                Y_array = reshape(Y_BC,1,num_points);
                Z_array = Rom.Physical_Displacement_Polynomial.evaluate_polynomial([X_array;Y_array]);
                num_dof = size(Z_array,1);
                Z_BC = reshape(Z_array,[num_dof,size(X_BC)]);
                Z_BC = permute(Z_BC,[2,3,1]);

                colour_data = ones(size(Z_BC,[1,2]));

                mesh(ax,Z_BC(:,:,plot_order(1)),Z_BC(:,:,plot_order(2)),Z_BC(:,:,plot_order(3)),colour_data,mesh_settings{:},"tag",manifold_name);
        end
        hold(ax,"off")
    end
%-------------------------------------------------------------------------%
    function ax = plot_orbit(ax,Dyn_Data,orbit_data,Plot_Settings,validation_manifold_plotting)
        LINE_WIDTH = 2;

        line_style = {"linewidth",LINE_WIDTH};
        orbit_style = "-";
        Rom = Dyn_Data.Dynamic_Model;

        
        solution_num = orbit_data(1);
        orbit_id = orbit_data(2);
        [orbit, validation_orbit] = Dyn_Data.get_orbit(solution_num,orbit_id,1);
        
        orbit_name = "o-" + join(string(Rom.Model.reduced_modes),",") + "-" + solution_num + "," + orbit_id;

        num_r_modes = length(Rom.Model.reduced_modes);
        if num_r_modes == 1
            orbit_style = ".-";
        end
        r_orbit = orbit.xbp(:,1:num_r_modes)';
        x_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_orbit);

        plot_order = Plot_Settings.coords;
        hold(ax,"on")
        plot3(ax,x_tilde(plot_order(1),:),x_tilde(plot_order(2),:),x_tilde(plot_order(3),:),orbit_style,line_style{:},"tag",orbit_name);
        hold(ax,"off")

        if isempty(validation_orbit)
            return
        end
        validation_orbit_name = "v" + orbit_name;
        Static_Data = load_static_data(Rom);
        Sol = Dyn_Data.load_solution(solution_num,"validation");
        Static_Data = Static_Data.add_validation_data(Sol.validation_modes);
        Rom = Reduced_System(Static_Data);

        dof = size(x_tilde,1);
        h_orbit = validation_orbit.h;
        num_points = size(h_orbit,2);
        x_hat = zeros(dof,num_points);
        Theta_Hat_Grad = Rom.Low_Frequency_Coupling_Gradient_Polynomial;
        for iPoint = 1:num_points
            x_hat_grad_orbit = Theta_Hat_Grad.evaluate_polynomial(r_orbit(:,iPoint));
            x_hat(:,iPoint) = x_tilde(:,iPoint) + x_hat_grad_orbit*h_orbit(:,iPoint);
        end

        hold(ax,"on")
        plot3(ax,x_hat(plot_order(1),:),x_hat(plot_order(2),:),x_hat(plot_order(3),:),"--",line_style{:},"tag",validation_orbit_name);
        hold(ax,"off")
        
        if validation_manifold_plotting
            % Static_Data = load_static_data(Rom);
            % Static_Data = Static_Data.add_validation_data(Sol.validation_modes,1);
            % Rom = Reduced_System(Static_Data);
            ax = plot_validation_manifold(ax,Rom,orbit_data,orbit.xbp',h_orbit,Plot_Settings);
        end
    end
%-------------------------------------------------------------------------%
    function ax = plot_validation_manifold(ax,Rom,orbit_data,r_orbit,h_orbit,Plot_Settings)
        
        NUM_DIMENSIONS = 2;
        H_LIM_SCALE_FACTOR = 1.5;
        PLOT_RESOLUTION = 101;
        MESH_ALPHA = 0.5;

        mesh_settings = {"EdgeColor","none","FaceColor","interp","FaceLighting","gouraud","FaceAlpha",MESH_ALPHA};
        manifold_name = "vm-" + join(string(Rom.Model.reduced_modes),",");

        plot_order = Plot_Settings.coords;
        %---------------------------------------
        r_transform = Rom.Model.reduced_eigenvectors'*Rom.Model.mass;

        switch NUM_DIMENSIONS
            case 1
                orbit_lim = [min(h_orbit(2,:)),max(h_orbit(2,:))].*H_LIM_SCALE_FACTOR;
                x = linspace(orbit_lim(1),orbit_lim(2),PLOT_RESOLUTION);
                
                
                num_t_points = size(r_orbit,2);
                hold(ax,"on")
                for iPoint = 1:num_t_points
                    r_i = r_orbit(:,iPoint);

                    theta_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_i);
                    displacement_gradient = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_i);
                    target_h = h_orbit(:,iPoint);
                    theta_hat = theta_tilde + displacement_gradient*target_h;

                    theta_plot = [theta_tilde,theta_hat];
                    theta_grad = diff(theta_plot,1,2);
                    theta_line = theta_tilde + theta_grad.*[-H_LIM_SCALE_FACTOR,H_LIM_SCALE_FACTOR];
                    
                    plot3(ax,theta_line(plot_order(1),:),theta_line(plot_order(2),:),theta_line(plot_order(3),:));
                end
                hold(ax,"off")
            case 2
                orbit_lim = [min(h_orbit,[],2),max(h_orbit,[],2)].*H_LIM_SCALE_FACTOR;

                num_points = PLOT_RESOLUTION^2;
                x = linspace(orbit_lim(1,1),orbit_lim(1,2),PLOT_RESOLUTION);
                y = linspace(orbit_lim(2,1),orbit_lim(2,2),PLOT_RESOLUTION);
                [X,Y] = meshgrid(x,y);

                X_array = reshape(X,1,num_points);
                Y_array = reshape(Y,1,num_points);
                h_array = [X_array;Y_array];


                num_t_points = size(r_orbit,2);

                Validation_Input = Rom.get_solver_inputs("h_prediction");

                num_r_modes = size(r_orbit,1)/2;
                r = r_orbit(1:num_r_modes,:);
                r_dot = r_orbit((num_r_modes+1):end,:);

                Eom_Input = Rom.get_solver_inputs("coco_backbone");
                

                
                hold(ax,"on")
                for iPoint = 1:num_t_points
                    r_i = r(:,iPoint);
                    r_dot_i = r_dot(:,iPoint);
                    
                    z_ddot = coco_eom(0,r_orbit(:,iPoint),0,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
                    r_ddot_i = z_ddot((num_r_modes+1):end,:);
                    [h_inertia,h_conv,h_stiff,h_force] = get_h_error_terms(r_i,r_dot_i,r_ddot_i,Validation_Input);
                    
                    
                    h = h_stiff\h_force;

                    displacement_gradient = Rom.Low_Frequency_Coupling_Gradient_Polynomial.evaluate_polynomial(r_i);
                    theta_tilde = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_i);
                    Z_array = theta_tilde + displacement_gradient*h;
                    plot3(Z_array(plot_order(1),:),Z_array(plot_order(2),:),Z_array(plot_order(3),:),"x")
                    % num_dof = size(Z_array,1);
                    % Z = reshape(Z_array,[num_dof,size(X)]);
                    % Z = permute(Z,[2,3,1]);
                    % colour_data = ones(size(Z,[1,2]));
                    % mesh(ax,Z(:,:,plot_order(1)),Z(:,:,plot_order(2)),Z(:,:,plot_order(3)),colour_data,mesh_settings{:},"tag",manifold_name);
                end
                hold(ax,"off")
        end

    end
end

