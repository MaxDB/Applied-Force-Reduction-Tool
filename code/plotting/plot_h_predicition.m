function ax = plot_h_predicition(Dyn_Data,type,solution_num,ax)
PLOT_BACKBONE = 1;

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
        energy_tilde = Dyn_Data.energy{1,solution_num};
        energy_hat = Dyn_Data.h_energy{1,solution_num};
        
        if ~exist("ax","var")
            figure
            ax = axes(gcf);
        end

        box(ax,"on")
        hold(ax,"on")
        


        p1 = plot(ax,frequency,energy_hat,'r-');
        if PLOT_BACKBONE
            p2 = plot(ax,frequency,energy_tilde,'k-');
            p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
        end
        hold(ax,"off")

        max_e_hat = max(energy_hat);
        max_e_tilde = max(energy_tilde);

        y_lim = min(max_e_tilde*3,max_e_hat);
        ylim(ax,[0,y_lim])

        p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
        

        xlabel(ax,"Frequency (rad/s)")
        ylabel(ax,"Energy")

    case "amplitude"
        r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        L_modes = Dyn_Data.validation_modes{1,solution_num};
        h_modes = [r_modes,L_modes];
        num_h_modes = length(h_modes);
        
        g_amp = Dyn_Data.corrected_low_modal_amplitude{1,solution_num};
        q_amp = Dyn_Data.low_modal_amplitude{1,solution_num};
        
        if ~exist("ax","var")
            figure
            tiledlayout("flow")
            for iMode = 1:num_h_modes
                ax{iMode,1} = nexttile;
                ax{iMode,2} = h_modes(iMode);
            end
            plotted_modes = h_modes;
        else 
            num_axes = size(ax,1);
            %check if all required r_modes are present]
            for iAx = 1:num_axes
                plotted_modes(iAx,1) = ax{iAx,2}; %#ok<AGROW>
            end

            neglected_modes = setdiff(h_modes,plotted_modes);
            num_neglected_modes = length(neglected_modes);
            for iMode = 1:num_neglected_modes
                ax{num_axes+iMode,1} = nexttile;
                ax{num_axes+iMode,2} = neglected_modes(iMode);
                plotted_modes(num_axes+iMode) = neglected_modes(iMode);
            end
        end

        for iMode = 1:num_h_modes
            mode = h_modes(iMode);
            ax_id = plotted_modes == mode;
            iAx = ax{ax_id,1};
            
            box(iAx,"on")
            
            hold(iAx,"on")
            p1 = plot(iAx,frequency,q_amp(iMode,:),'r-');
            p2 = plot(iAx,frequency,g_amp(iMode,:),'k-');
            hold(iAx,"off")

            max_q = max(q_amp(iMode,:));
            max_g = max(g_amp(iMode,:));
            
            y_lim = min(max_g*3,max_q);
            ylim(iAx,[0,y_lim])

            p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
            p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

            xlabel(iAx,"Frequency (rad/s)")
            ylabel(iAx,"Q_{" + h_modes(iMode) + "}")
        end
end
end