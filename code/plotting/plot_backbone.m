function ax = plot_backbone(Dyn_Data,type,solution_num,ax)
if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end

if isstring(solution_num)
    if solution_num == "last"
        solution_num = Dyn_Data.num_solutions;
    end
end

frequency = Dyn_Data.frequency{1,solution_num};
num_orbits = size(frequency,2);
orbit_labels = 1:num_orbits;
orbit_ids = "(" + solution_num + "," + orbit_labels + ")";
data_tip_row = dataTipTextRow("ID",orbit_ids);

switch type
    case "energy"
        energy = Dyn_Data.energy{1,solution_num};

        if ~exist("ax","var")
            figure
            ax = axes(gcf);
        end

        box(ax,"on")

        hold(ax,"on")
        p = plot(ax,frequency,energy);
        hold(ax,"off")

        
        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

        xlabel(ax,"Frequency (rad/s)")
        ylabel(ax,"Energy")
    case "amplitude"
        amplitude = Dyn_Data.amplitude{1,solution_num};
        solution_type = Dyn_Data.solution_types{1,solution_num};
        num_modes = size(amplitude,1);

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
            p = plot(ax{ax_id,1},frequency,amplitude(iMode,:));
            hold(ax{ax_id,1},"off")

            xlabel(ax{ax_id,1},"Frequency (rad/s)")
            ylabel(ax{ax_id,1},"R_{" + r_modes(iMode) + "}")

            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
        end
end
end