function ax = plot_orbit(Dyn_Data,type,solution_num,orbit_num,varargin)

LINE_WIDTH = 1;
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
colour_num = 1;
normalise = 1;
orbit_shift = 1;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case {"colour","color"}
            colour_num = keyword_values{arg_counter};
        case {"normalise"}
            normalise = keyword_values{arg_counter};
        case {"shift_orbit"}
            orbit_shift = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",line_colour};


if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
orbit = Dyn_Data.get_orbit(solution_num,orbit_num);
time = orbit.tbp';

if normalise
    period = orbit.T;
    time = time/period;
end

switch type
    case "displacement"
        displacement = orbit.xbp';
        num_modes = size(displacement,1)/2;
        displacement = displacement(1:num_modes,:);

        Solution_Type = Dyn_Data.solution_types{1,solution_num};
  
         

        if Solution_Type.model_type == "fom"
            r_modes = 1:num_modes;
        else
            r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        end


        if isempty(ax)
            figure
            tiledlayout("flow")
            for iMode = 1:num_modes
                ax{iMode,1} = nexttile; %#ok<*AGROW>
                ax{iMode,2} = r_modes(iMode);
            end
            plotted_modes = r_modes;
        else
            num_axes = size(ax,1);
            %check if all required r_modes are present]
            for iAx = 1:num_axes
                plotted_modes(iAx,1) = ax{iAx,2}; 
            end

            neglected_modes = setdiff(r_modes,plotted_modes);
            num_neglected_modes = length(neglected_modes);
            for iMode = 1:num_neglected_modes
                ax{num_axes+iMode,1} = nexttile;
                ax{num_axes+iMode,2} = neglected_modes(iMode);
                plotted_modes(num_axes+iMode) = neglected_modes(iMode);
            end
        end



        for iMode = 1:num_modes
            mode = r_modes(iMode);
            ax_id = plotted_modes == mode;
            box(ax{ax_id,1},"on")

            hold(ax{ax_id,1},"on")
            plot(ax{ax_id,1},time,displacement(iMode,:),line_plot_settings{:});
            hold(ax{ax_id,1},"off")

            xlabel(ax{ax_id,1},"Time (s)")
            ylabel(ax{ax_id,1},"r_{" + r_modes(iMode) + "}")

        end
    case "physical displacement"
        displacement = orbit.xbp';
        num_modes = size(displacement,1)/2;
        displacement = displacement(1:num_modes,:);

        Rom = Dyn_Data.Dynamic_Model;
        node_map = Rom.Model.node_mapping;
        dof_bc = Dyn_Data.Additional_Output.dof;
        dof = node_map(node_map(:,1) == dof_bc,2);

        x_phy = Rom.expand(displacement);
        x_dof = x_phy(dof,:);
        
        if isempty(ax)
            figure;
            ax = gca;
            box(ax,"on")
            ylabel(ax,"x_{" + dof_bc + "} (m)")
            if normalise
                xlabel(ax,"Normalised Time")
                xlim(ax,[0,1])
            else
                xlabel(ax,"Time (s)")
            end
        end
        
        if orbit_shift
            [time,x_dof] = shift_orbit(time,x_dof);
        end

        hold(ax,"on")
        plot(time,x_dof,line_plot_settings{:})
        hold(ax,"off")

end
end