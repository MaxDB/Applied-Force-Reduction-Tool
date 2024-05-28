function ax = plot_orbit(Dyn_Data,type,solution_num,orbit_num,ax)
if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
orbit = Dyn_Data.get_orbit(solution_num,orbit_num);
time = orbit.tbp';

switch type
    case "displacement"
        displacement = orbit.xbp';
        num_modes = size(displacement,1)/2;
        displacement = displacement(1:num_modes,:);

        solution_type = Dyn_Data.solution_types{1,solution_num};
  
         

        if solution_type == "fom: backbone"
            r_modes = 1:num_modes;
        else
            r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        end


        if ~exist("ax","var")
            figure
            tiledlayout("flow")
            for iMode = 1:num_modes
                ax{iMode,1} = nexttile;
                ax{iMode,2} = r_modes(iMode);
            end
            plotted_modes = r_modes;
        else
            num_axes = size(ax,1);
            %check if all required r_modes are present]
            for iAx = 1:num_axes
                plotted_modes(iAx,1) = ax{iAx,2}; %#ok<AGROW>
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
            plot(ax{ax_id,1},time,displacement(iMode,:));
            hold(ax{ax_id,1},"off")

            xlabel(ax{ax_id,1},"Time (s)")
            ylabel(ax{ax_id,1},"r_{" + r_modes(iMode) + "}")

        end
end
end