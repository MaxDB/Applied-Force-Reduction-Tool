function ax = plot_backbone(Dyn_Data,type,solution_num,ax)
PLOT_BIFURCATIONS = 1;
PLOT_PERIODICITY = 0;
LINE_STYLE = [":","-"]; %[unstable,stable]
LINE_WIDTH = 1;
LINE_COLOUR = "k";
BIFURCATION_MARKER = ["o","o","^","x"];
BIFURCATION_TYPE = ["BP","PD","NS","SN"];
BIFURCATION_SIZE = [6,6,6,0];
STABILITY_LIMIT = 1.005;

%-------------------------------------------------------------------------%
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

line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",LINE_COLOUR};

if PLOT_BIFURCATIONS
    bifurcations = Dyn_Data.bifurcations{1,solution_num};
    bifurcation_types = fields(bifurcations);
    num_bifurcation_types = size(bifurcation_types,1);

    stability = Dyn_Data.stability{1,solution_num};
    [index_ranges,range_stability] = get_stability_boundaries(stability<STABILITY_LIMIT,bifurcations);
    num_sections = size(range_stability,2);

    bifurcation_plot_settings = cell(num_bifurcation_types,1);
    for iType = 1:num_bifurcation_types
        type_plot_settings = {"LineStyle","none","LineWidth",LINE_WIDTH,"Color",LINE_COLOUR,...
            "Marker",BIFURCATION_MARKER(iType),"MarkerSize",BIFURCATION_SIZE(iType)};
        switch iType
            case 1
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = "w";
            case 2
                type_plot_settings{end+1} = "MarkerFaceColor";
                type_plot_settings{end+1} = "w";
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = LINE_COLOUR; %#ok<*AGROW>
            case 3
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = "w";
        end
        bifurcation_plot_settings{iType,1} =  type_plot_settings;
    end
end

periodicity_error = Dyn_Data.periodicity_error{1,solution_num};
if isempty(periodicity_error)
    PLOT_PERIODICITY = 0;
end

switch type
    case {"energy","physical amplitude"}
        switch type
            case "energy"
                energy = Dyn_Data.energy{1,solution_num};
            case "physical amplitude"
                energy = Dyn_Data.additional_dynamic_output{1,solution_num};
        end
        

        if ~exist("ax","var")
            figure
            ax = axes(gcf);
        end

        box(ax,"on")

        hold(ax,"on")
        if PLOT_BIFURCATIONS
            for iSection = 1:num_sections
                stab = range_stability(iSection);
                index_range = index_ranges(iSection,1):index_ranges(iSection,2);
                p = plot(ax,frequency(index_range),energy(index_range),'LineStyle',LINE_STYLE(stab+1),line_plot_settings{:});

                data_tip_row = dataTipTextRow("ID",orbit_ids(index_range));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;


                data_tip_row = dataTipTextRow("Stability",stability(index_range));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                if PLOT_PERIODICITY
                    data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end

            end

            for iType = 1:num_bifurcation_types
                bifurcation_type = bifurcation_types{iType,1};
                bifurcation_index = bifurcations.(bifurcation_type);
                if isempty(bifurcation_index)
                    continue
                end
                if BIFURCATION_SIZE(iType) == 0
                    continue
                end

                type_plot_settings = bifurcation_plot_settings{iType,1};
                p = plot(ax,frequency(bifurcation_index),energy(bifurcation_index),type_plot_settings{:});

                data_tip_row = dataTipTextRow("ID",orbit_ids(bifurcation_index));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                data_tip_row = dataTipTextRow("Stability",stability(bifurcation_index));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                data_tip_row = dataTipTextRow("Type",repelem(BIFURCATION_TYPE(iType),size(bifurcation_index,1),size(bifurcation_index,2)));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                

                if PLOT_PERIODICITY
                    data_tip_row = dataTipTextRow("Periodicity",periodicity_error(bifurcation_index));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end
            end

        else
            p = plot(ax,frequency,energy,'LineStyle',LINE_STYLE(2),line_plot_settings{:});
            data_tip_row = dataTipTextRow("ID",orbit_ids);
            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

            if PLOT_PERIODICITY
                data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
            end
        end
        hold(ax,"off")



        xlabel(ax,"Frequency (rad/s)")
        switch type
            case "energy"
                ylabel(ax,"Energy")
            case "physical amplitude"
                output_dof = Dyn_Data.Additional_Output.dof;
                ylabel(ax,"X_{" + output_dof + "}");
        end
        
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
            if PLOT_BIFURCATIONS
                for iSection = 1:num_sections
                    stab = range_stability(iSection);
                    index_range = index_ranges(iSection,1):index_ranges(iSection,2);
                    p = plot(ax{ax_id,1},frequency(index_range),amplitude(iMode,index_range),'LineStyle',LINE_STYLE(stab+1),line_plot_settings{:});

                    data_tip_row = dataTipTextRow("ID",orbit_ids(index_range));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;


                    data_tip_row = dataTipTextRow("Stability",stability(index_range));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    if PLOT_PERIODICITY
                        data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end

                end

                for iType = 1:num_bifurcation_types
                    bifurcation_type = bifurcation_types{iType,1};
                    bifurcation_index = bifurcations.(bifurcation_type);
                    if isempty(bifurcation_index)
                        continue
                    end
                    if BIFURCATION_SIZE(iType) == 0
                        continue
                    end

                    type_plot_settings = bifurcation_plot_settings{iType,1};
                    p = plot(ax{ax_id,1},frequency(bifurcation_index),amplitude(iMode,bifurcation_index),type_plot_settings{:});

                    data_tip_row = dataTipTextRow("ID",orbit_ids(bifurcation_index));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    data_tip_row = dataTipTextRow("Stability",stability(bifurcation_index));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    data_tip_row = dataTipTextRow("Type",repelem(BIFURCATION_TYPE(iType),size(bifurcation_index,1),size(bifurcation_index,2)));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    if PLOT_PERIODICITY
                        data_tip_row = dataTipTextRow("Periodicity",periodicity_error(bifurcation_index));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end
                end
            else
                p = plot(ax{ax_id,1},frequency,amplitude(iMode,:),'LineStyle',LINE_STYLE(2),line_plot_settings{:});
                data_tip_row = dataTipTextRow("ID",orbit_ids);
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                if PLOT_PERIODICITY
                    data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end
            end
            hold(ax{ax_id,1},"off")

            xlabel(ax{ax_id,1},"Frequency (rad/s)")
            ylabel(ax{ax_id,1},"R_{" + r_modes(iMode) + "}")


        end
end
end