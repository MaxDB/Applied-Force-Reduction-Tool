function ax = plot_floquet_multipliers(Dyn_Data,solution_num,orbit_num,varargin)
LINE_WIDTH = 2;

%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
colour_num = 1;
tag = "";

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case {"colour","color"}
            colour_num = keyword_values{arg_counter};
        case "tag"
            tag = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"MarkerSize",10,"Color",line_colour,"DisplayName","r","Tag",tag};

if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
Orbit = Dyn_Data.get_orbit(solution_num,orbit_num);
Rom = Dyn_Data.Dynamic_Model;
Model = Rom.Model;

Sol = Dyn_Data.load_solution(solution_num);
Sol_Type = Dyn_Data.solution_types{1,solution_num};


floquet_multipliers = Orbit.po_test.la;


if isempty(ax)
    fig = figure;
    ax = axes;

    theta = linspace(-pi,pi,50);
    x = cos(theta);
    y = sin(theta);
    hold(ax,"on")
    plot(ax,x,y)
    hold(ax,"off")

    ax.DataAspectRatio = [1,1,1];
    box(ax,"on")
end



mult_x = real(floquet_multipliers);
mult_y = imag(floquet_multipliers);

hold(ax,"on")
plot(ax,mult_x,mult_y,"x",line_plot_settings{:})
hold(ax,"off")



end