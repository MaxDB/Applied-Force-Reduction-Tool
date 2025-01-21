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
    if class(keyword_args{arg_counter}) == "Dynamic_Dataset"
        system_counter = system_counter + 1;
        dyn_data_names{1,system_counter} = keyword_args{arg_counter};
        solution_index{1,system_counter} = keyword_values{arg_counter};
        continue
    end
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

if num_solutions == 1
    colour_numbers = 0;
else
    colour_numbers = 1:8;
end
ax = [];
legend_lines = zeros(1,num_solutions);
legend_modes = cell(1,num_solutions);
for iSol = 1:num_solutions
    colour_sol = colour_numbers(iSol);
    data_name = dyn_data_names{1,iSol};
    Dyn_Data = initalise_dynamic_data(data_name);
    

    system_sols = solution_index{1,iSol};
    if isstring(system_sols) && system_sols == "all"
        system_sols = 1:size(Dyn_Data,1);
    end
    num_system_sols = length(system_sols);
    for jSol = 1:num_system_sols
        Solution = Dyn_Data.load_solution(system_sols(jSol));
        Sol_Type = Solution.Solution_Type;


        switch Sol_Type.orbit_type
            case "forced"
                
                if isfield(Sol_Type,"amplitude")
                    amplitude = Sol_Type.amplitude;
                    if ~exist("amplitude_map","var")
                        amplitude_map = dictionary;
                        amplitude_map(round(amplitude)) = 1;
                    else
                       if ~isKey(amplitude_map,round(amplitude))
                            amplitude_map(round(amplitude)) = numEntries(amplitude_map) + 1;
                       end
                    end
                    colour_sol = amplitude_map(round(amplitude));
                else
                    colour_sol = "grey";
                end
            case "free"
                if num_solutions == 1
                    colour_sol = 0;
                end
        end
        
        if Sol_Type.validated == 0
            validation(iSol) = 0;
        end

        if validation(iSol)
            if colour_sol == 0
                h_colour_number = 1;
            else
                h_colour_number = colour_sol;
            end

            ax = plot_h_predicition(Dyn_Data,type,system_sols(jSol),"axes",ax,"colour",h_colour_number);
        else
            ax = plot_backbone(Dyn_Data,type,system_sols(jSol),"axes",ax,"colour",colour_sol);
            if Sol_Type.orbit_type == "free"
                if iscell(ax)
                    num_axes = size(ax,1);
                    for iAx = 1:num_axes
                        uistack(ax{iAx,1}.Children(1),"bottom")
                    end
                else
                    uistack(ax(1).Children(1),"bottom")
                end

            end
        end
    end
    if PLOT_LEGEND
        if iscell(ax)
            ax_legend = ax{1,1};
        else
            ax_legend = ax;
        end
        hold(ax_legend,"on")
        line = ax_legend.Children(1);
        point = [line.XData(1),line.YData(1)];
        line_width = line.LineWidth;
        line_colour = get_plot_colours(colour_numbers(iSol));
        p = plot(ax_legend,point(1),point(2),"-","LineWidth",line_width,"Color",line_colour);
        hold(ax_legend,"off")
        legend_lines(iSol) = p;
        
        switch Sol_Type.model_type
            case "rom"
                reduced_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
            case "fom"
                reduced_modes = 1:Dyn_Data.Dynamic_Model.Model.num_dof;
        end
        legend_modes{1,iSol} = reduced_modes;
    end
end

if PLOT_LEGEND
    
    legend_names = cellfun(@(x) "\{" + join(string(x),", ") + "\}",legend_modes);
    legend(legend_lines,legend_names,"Location","best")
end
end


