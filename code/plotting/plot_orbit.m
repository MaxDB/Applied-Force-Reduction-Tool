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
validated = 0;

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
        case {"validated"}
            validated = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",line_colour};

if validated
    validated_line_colour = get_plot_colours(colour_num+1);
    validated_line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",validated_line_colour};
end


if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
[orbit,validated_orbit] = Dyn_Data.get_orbit(solution_num,orbit_num,validated);
time = orbit.tbp';


if normalise
    period = orbit.T;
    time = time/period;
end

switch type
    case {"displacement","phase"}
        displacement = orbit.xbp';
        num_r_modes = size(displacement,1)/2;
        
        displacement = displacement(1:num_r_modes,:);

        Solution_Type = Dyn_Data.solution_types{1,solution_num};
  
         

        if Solution_Type.model_type == "fom"
            r_modes = 1:num_r_modes;
        else
            r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        end
        
        if validated
            L_modes = Dyn_Data.validation_modes{1,solution_num};

            physical_displacement = Dyn_Data.Dynamic_Model.expand(displacement);
            Model = Dyn_Data.Dynamic_Model.Model;
            mass = Model.mass;
            stiffness = Model.stiffness;

            [full_evecs,~] = eigs(stiffness,mass,max(L_modes),"smallestabs");
            % full_evals = diag(full_evals);
            new_L_modes = 1:max(L_modes);
            new_L_modes(ismember(new_L_modes,r_modes)) = [];

            L_evecs = full_evecs(:,new_L_modes);
            % r_evecs = Model.reduced_eigenvectors;

            g = L_evecs'*mass*physical_displacement;
            
            displacement = [displacement;g];

            h_displacement = validated_orbit.h;
            h_displacement = displacement + h_displacement;


            h_modes = [r_modes,L_modes];
            plot_modes = h_modes;
        else
            plot_modes = r_modes;
        end

        num_modes = size(plot_modes,2);

        switch type
            case "displacement"

                if isempty(ax)
                    figure
                    tiledlayout("flow")
                    for iMode = 1:num_modes
                        ax{iMode,1} = nexttile; %#ok<*AGROW>
                        ax{iMode,2} = plot_modes(iMode);
                    end
                    plotted_modes = plot_modes;
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
                    mode = plot_modes(iMode);
                    ax_id = plotted_modes == mode;
                    box(ax{ax_id,1},"on")
                    

                    hold(ax{ax_id,1},"on")
                    plot(ax{ax_id,1},time,displacement(iMode,:),line_plot_settings{:});
                    if validated
                        plot(ax{ax_id,1},time,h_displacement(iMode,:),validated_line_plot_settings{:});
                    end
                    hold(ax{ax_id,1},"off")

                    xlabel(ax{ax_id,1},"Time (s)")
                    ylabel(ax{ax_id,1},"r_{" + plot_modes(iMode) + "}")
                end
                title(ax{1,1},"(" + solution_num + "," + orbit_num + ")")
            case "phase"
                if isempty(ax)
                    figure;
                    ax = gca;
                    box(ax,"on")
                    title(ax,"(" + solution_num + "," + orbit_num + ")")
                    
                    if validated
                        xlabel("r") 
                        ylabel("s") 
                    else
                        xlabel("r") 
                        ylabel("r_dot") 
                    end
                end

                if validated
                    hold(ax,"on")
                    plot(ax,displacement(1,:),displacement(2,:),line_plot_settings{:})
                    plot(ax,h_displacement(1,:),h_displacement(2,:),validated_line_plot_settings{:})
                    hold(ax,"off")
                    xlabel("r")
                    ylabel("s")
                else
                    plot(ax,displacement(1,:),displacement(2,:),line_plot_settings{:})
                end

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