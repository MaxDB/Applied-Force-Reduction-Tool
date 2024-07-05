function ax = plot_h_predicition(Dyn_Data,type,solution_num,varargin)
PLOT_STABILITY = 1;
STABILITY_LIMIT = 1.005;

LINE_STYLE = [":","-"]; %[unstable,stable]
LINE_WIDTH = 1.5;


%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
colour_num = 1;
add_backbone = 1;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case {"colour","color"}
            colour_num = keyword_values{arg_counter};
        case {"backbone"}
            add_backbone = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end

if isstring(solution_num)
    if solution_num == "last"
        solution_num = Dyn_Data.num_solutions;
    end
end
%-------------------------------------------------------------------------%

frequency = Dyn_Data.frequency{1,solution_num};
num_orbits = size(frequency,2);


orbit_labels = 1:num_orbits;
orbit_ids = "(" + solution_num + "," + orbit_labels + ")";

validation_modes = Dyn_Data.validation_modes{1,solution_num};
mode_string = "[" + join(string(validation_modes),", ") + "]";

line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",line_colour};

if PLOT_STABILITY
    stability = Dyn_Data.h_stability{1,solution_num}; %#ok<*UNRCH>
    [index_ranges,range_stability] = get_stability_boundaries(stability<STABILITY_LIMIT);
    num_sections = size(range_stability,2);

    bb_stability = Dyn_Data.stability{1,solution_num};
    [bb_index_ranges,bb_range_stability] = get_stability_boundaries(bb_stability<STABILITY_LIMIT);
    bb_num_sections = size(bb_range_stability,2);

end

switch type
    case "energy"
        energy_tilde = Dyn_Data.energy{1,solution_num};
        energy_hat = Dyn_Data.h_energy{1,solution_num};
        
        if isempty(ax)
            figure
            ax = axes(gcf);
        end
        
        
        
        box(ax,"on")
        hold(ax,"on")

        % if add_backbone
        %     ax = plot_backbone(Dyn_Data,type,solution_num,"axes",ax,"colour",0);
        % end
        
        if PLOT_STABILITY
            for iSection = 1:num_sections
                stab = range_stability(iSection);
                index_range = index_ranges(iSection,1):index_ranges(iSection,2);
                p1 = plot(ax,frequency(index_range),energy_hat(index_range),'LineStyle',LINE_STYLE(stab+1),line_plot_settings{:});

                data_tip_row_id = dataTipTextRow("ID",orbit_ids(index_range));
                p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

                data_tip_row_modes = dataTipTextRow("Modes",repelem(mode_string,1,num_orbits));
                p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_modes;

                data_tip_row_stab = dataTipTextRow("Stability",stability(index_range));
                p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_stab;
            end
        else
            p1 = plot(ax,frequency,energy_hat,'LineStyle',LINE_STYLE(2),line_plot_settings{:});
            data_tip_row = dataTipTextRow("ID",orbit_ids);
            p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

            data_tip_row_modes = dataTipTextRow("Modes",repelem(mode_string,1,num_orbits));
            p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_modes;
        end

       

        
        hold(ax,"off")

        max_e_hat = max(energy_hat);
        max_e_tilde = max(energy_tilde);

        y_lim = min(max_e_tilde*3,max_e_hat);
        ylim(ax,[0,y_lim])

        xlabel(ax,"Frequency (rad/s)")
        ylabel(ax,"Energy")

    case "amplitude"
        r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        L_modes = Dyn_Data.validation_modes{1,solution_num};
        h_modes = [r_modes,L_modes];

        num_h_modes = length(h_modes);
        
        g_amp = Dyn_Data.corrected_low_modal_amplitude{1,solution_num};
        q_amp = Dyn_Data.low_modal_amplitude{1,solution_num};
        
        if isempty(ax)
            figure
            tiledlayout("flow")
            for iMode = 1:num_h_modes
                ax{iMode,1} = nexttile; %#ok<AGROW>
                ax{iMode,2} = h_modes(iMode); %#ok<AGROW>
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
                ax{num_axes+iMode,1} = nexttile; %#ok<AGROW>
                ax{num_axes+iMode,2} = neglected_modes(iMode); %#ok<AGROW>
                plotted_modes(num_axes+iMode) = neglected_modes(iMode);
            end
        end

        for iMode = 1:num_h_modes
            mode = h_modes(iMode);
            ax_id = plotted_modes == mode;
            iAx = ax{ax_id,1};
            
            box(iAx,"on")

            hold(iAx,"on")
            

            if PLOT_STABILITY
                for iSection = 1:num_sections
                    stab = range_stability(iSection);
                    index_range = index_ranges(iSection,1):index_ranges(iSection,2);
                    p1 = plot(iAx,frequency(index_range),q_amp(iMode,index_range),'LineStyle',LINE_STYLE(stab+1),line_plot_settings{:});

                    data_tip_row_id = dataTipTextRow("ID",orbit_ids(index_range));
                    p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

                    data_tip_row_modes = dataTipTextRow("Modes",repelem(mode_string,1,num_orbits));
                    p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_modes;

                    data_tip_row_stab = dataTipTextRow("Stability",stability(index_range));
                    p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_stab;
                end
            else
                p1 = plot(iAx,frequency,q_amp(iMode,:),line_plot_settings{:}); 

                data_tip_row_id = dataTipTextRow("ID",orbit_ids);
                p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

                data_tip_row_modes = dataTipTextRow("Modes",repelem(mode_string,1,num_orbits));
                p1.DataTipTemplate.DataTipRows(end+1) = data_tip_row_modes;
            end

            if add_backbone
                if PLOT_STABILITY
                    for iSection = 1:bb_num_sections
                        stab = bb_range_stability(iSection);
                        index_range = bb_index_ranges(iSection,1):bb_index_ranges(iSection,2);
                        p2 = plot(iAx,frequency(index_range),g_amp(iMode,index_range),"color",[0,0,0],'LineStyle',LINE_STYLE(stab+1),line_plot_settings{1:2});

                        data_tip_row_id = dataTipTextRow("ID",orbit_ids(index_range));
                        p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

                        data_tip_row_modes = dataTipTextRow("Modes",repelem(mode_string,1,num_orbits));
                        p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row_modes;

                        data_tip_row_stab = dataTipTextRow("Stability",stability(index_range));
                        p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row_stab;
                    end
                else
                    p2 = plot(iAx,frequency,g_amp(iMode,:),'k-',line_plot_settings{1:2}); 

                    data_tip_row_id = dataTipTextRow("ID",orbit_ids);
                    p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

                    data_tip_row_modes = dataTipTextRow("Modes",repelem(mode_string,1,num_orbits));
                    p2.DataTipTemplate.DataTipRows(end+1) = data_tip_row_modes;
                end
            end

            hold(iAx,"off")

            max_q = max(q_amp(iMode,:));
            max_g = max(g_amp(iMode,:));
            
            y_lim = min(max_g*3,max_q);
            if y_lim ~= 0
                ylim(iAx,[0,y_lim])
            end

            xlabel(iAx,"Frequency (rad/s)")
            ylabel(iAx,"Q_{" + h_modes(iMode) + "}")
        end

        
end
end