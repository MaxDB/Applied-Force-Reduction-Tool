function compare_solutions(type,varargin)
PLOT_LEGEND = 1;

%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

dyn_data_names = cell(0,1);
solution_index = cell(0,1);
system_counter = 0;
for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "validation"
            validation = keyword_values{arg_counter};
        otherwise
            system_counter = system_counter + 1;
            dyn_data_names{1,system_counter} = keyword_args{arg_counter};
            solution_index{1,system_counter} = keyword_values{arg_counter};
            
    end
end
%-------------------------------------------------------------------------%
num_solutions = system_counter;
if ~exist("validation","var")
    validation = zeros(1,num_solutions);
elseif isscalar(validation)
    validation = validation * ones(1,num_solutions);
end

colour_numbers = 1:num_solutions;

ax = [];
legend_lines = zeros(1,num_solutions);
legend_modes = cell(1,num_solutions);
for iSol = 1:num_solutions
    data_name = dyn_data_names{1,iSol};
    Dyn_Data = initalise_dynamic_data(data_name);

    system_sols = solution_index{1,iSol};
    if isstring(system_sols) && system_sols == "all"
        system_sols = 1:size(Dyn_Data,1);
    end
    num_system_sols = length(system_sols);
    for jSol = 1:num_system_sols
        if validation(num_solutions)
            ax = plot_h_predicition(Dyn_Data,type,system_sols(jSol),"axes",ax,"colour",colour_numbers(iSol));
        else
            ax = plot_backbone(Dyn_Data,type,system_sols(jSol),"axes",ax,"colour",colour_numbers(iSol));
        end
    end
    if PLOT_LEGEND
        hold(ax,"on")
        line = ax.Children(1);
        point = [line.XData(1),line.YData(1)];
        line_width = line.LineWidth;
        line_colour = get_plot_colours(colour_numbers(iSol));
        p = plot(ax,point(1),point(2),"-","LineWidth",line_width,"Color",line_colour);
        hold(ax,"off")
        legend_lines(iSol) = p;
        reduced_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        legend_modes{1,iSol} = reduced_modes;
    end
end

if PLOT_LEGEND
    
    legend_names = cellfun(@(x) "\{" + join(string(x),", ") + "\}",legend_modes);
    legend(legend_lines,legend_names,"Location","best")
end
end


