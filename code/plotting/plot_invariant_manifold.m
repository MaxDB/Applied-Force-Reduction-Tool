function ax = plot_invariant_manifold(Dyn_Data,sol_num,Plotting_Opts,varargin)
%-------------------------------------------------------------------------%
points_per_orbit = 500;
plotting_coords = [1,3,2];
orbit_spacing = 0.002;
% add this into plotting options
%-------------------------------------------------------------------------%
mesh_alpha = 1;
mesh_settings = {"EdgeColor","none","LineWidth",0.01,"FaceColor",get_plot_colours(5),"FaceLighting","gouraud","FaceAlpha",mesh_alpha,"Tag","invariant_manifold"};
grid_line_style = {"Color",[0.25,0.25,0.25,0.2],"LineWidth",0.2,"Tag","grid_line"};
outline_style = {"Color",[0,0,0],"LineWidth",1,"Tag","outline"};
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
limit_type = {};
limit_value = {};
limit_counter = 0;

for arg_counter = 1:num_args/2
    keyword_arg = keyword_args{arg_counter};
    switch keyword_arg
        case "axes"
            ax = keyword_values{arg_counter};
        otherwise
            if endsWith(keyword_arg,"limit")
                limit_counter = limit_counter + 1;
                arg_array = split(keyword_arg," ");
                limit_type{limit_counter} = arg_array(1);
                limit_value{limit_counter} = keyword_values{arg_counter}; %#ok<*AGROW>
                continue
            end
            erorr("Invalid option: '" + keyword_args{arg_counter} + "'")
    end
end
%-------------------------------------------------------------------------%
if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
Rom = Dyn_Data.Dynamic_Model;
Model = Rom.Model;

evec = Model.reduced_eigenvectors;
num_modes = size(evec,2);
num_dof = Model.num_dof;

num_coords = size(plotting_coords,2);
disp_coord_index = plotting_coords <= num_dof;
disp_coords = plotting_coords(disp_coord_index);
vel_coords = plotting_coords(~disp_coord_index) - num_dof;

Sol = Dyn_Data.load_solution(sol_num);
num_orbits = Sol.num_orbits;

sol_points = zeros(num_orbits,points_per_orbit+1,3);
orbit_counter = 0;
for iOrbit = 1:num_orbits
    orbit = Dyn_Data.get_orbit(sol_num,iOrbit);

    period = orbit.T;
    freq = 2*pi/period;

    z = orbit.xbp';
    r = z(1:num_modes,:);    

    limit_reached = 0;
    for iLim = 1:limit_counter
        switch limit_type{iLim}
            case "energy"
                E = max(Rom.Potential_Polynomial.evaluate_polynomial(r));
                if E > limit_value{iLim}
                    limit_reached = 1;
                    continue
                end
            case "frequency"
                freq_range = limit_value{iLim};
                if freq < freq_range(1) || freq > freq_range(2)
                    limit_reached = 1;
                    continue
                end
            case "id"
                id_range = limit_value{iLim};
                if isscalar(id_range)
                    id_range = [1,id_range];
                end
                if iOrbit < id_range(1) || iOrbit > id_range(2)
                    limit_reached = 1;
                    continue
                end

            case "amplitude"

            case "disp"
                if max(abs(x(1,:))) > max_x1
                    limit_reached = 1;
                    continue
                end
            otherwise
                error("unknown limit type: '" + limit_type{iLim} + "'")
        end
    end
        
    if limit_reached
        continue
    end

    orbit_counter = orbit_counter + 1;

    t = orbit.tbp';
    num_t_points = size(t,2);
    t_lin = linspace(0,period,num_t_points);
    z_lin = zeros(2*num_modes,num_t_points);
    for iState = 1:(2*num_modes)
        z_lin(iState,:) = interp1(t,z(iState,:),t_lin);
    end
    z_interp = interpft(z_lin,points_per_orbit,2);
    z_interp = [z_interp,z_interp(:,1)]; 

    r_interp = z_interp(1:num_modes,:);
    r_dot_interp = z_interp((1:num_modes)+num_modes,:);

    x = Rom.Physical_Displacement_Polynomial.evaluate_polynomial(r_interp,disp_coords);
    x_dot = Rom.get_physical_velocity(r_interp,r_dot_interp,vel_coords);

    x_all = zeros(num_coords,points_per_orbit + 1);
    x_all(disp_coord_index,:) = x;
    x_all(~disp_coord_index,:) = x_dot;

    sol_points(orbit_counter,:,:) = x_all';
end
sol_points((orbit_counter+1):end,:,:) = [];


if isempty(ax)
    fig = figure;
    ax = axes(fig);
    box(ax,"on")
end

hold(ax,"on")
num_orbits = size(sol_points,1);
if num_orbits < 2
    return
end
mesh(ax,sol_points(:,:,1),sol_points(:,:,2),sol_points(:,:,3),mesh_settings{:})
plot3(ax,sol_points(1,:,1),sol_points(1,:,2),sol_points(1,:,3),outline_style{:})
last_orbit = zeros(3,points_per_orbit+1);
for iOrbit = 1:(num_orbits-1)
    orbit = squeeze(sol_points(iOrbit,:,:))';
    [~,max_distance] = orbit_distance(orbit,last_orbit);
    if max_distance < orbit_spacing
        continue
    end
    plot3(ax,sol_points(iOrbit,:,1),sol_points(iOrbit,:,2),sol_points(iOrbit,:,3),grid_line_style{:})
    last_orbit = orbit;
end
plot3(ax,sol_points(end,:,1),sol_points(end,:,2),sol_points(end,:,3),outline_style{:})
hold(ax,"off")

coord_names = strings(num_coords,1);
for iCoord = 1:num_coords
    if disp_coord_index(iCoord)
        coord_names(iCoord) = "$x_{" + plotting_coords(iCoord) + "}$";
    else
        coord_names(iCoord) = "$\dot{x}_{" + (plotting_coords(iCoord) - num_dof) + "}$";
    end
end
xlabel(ax,coord_names(1),"Interpreter","latex")
ylabel(ax,coord_names(2),"Interpreter","latex")
zlabel(ax,coord_names(3),"Interpreter","latex")
end

function [min_distance,max_distance] = orbit_distance(orbit_one,orbit_two)
    displacement = orbit_two - orbit_one;
    distance = sum(displacement.^2,1);
    min_distance = min(distance);
    max_distance = max(distance);
end