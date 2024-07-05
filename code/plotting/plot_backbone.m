function ax = plot_backbone(Dyn_Data,type,solution_num,varargin)
PLOT_BIFURCATIONS = 1;
PLOT_SPECIAL_POINT = 0;
PLOT_BB_SN = 0;

PLOT_PERIODICITY = 1;
STABILITY_LIMIT = 1.005;

LINE_STYLE = [":","-"]; %[unstable,stable]
LINE_WIDTH = 1.5;
BIFURCATION_MARKER = ["o","o","^","x"];
BIFURCATION_TYPE = ["BP","PD","NS","SN"];
BIFURCATION_SIZE = [6,0,6,6];


%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
colour_num = 1;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case {"colour","color"}
            colour_num = keyword_values{arg_counter};
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

line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",line_colour};

if PLOT_BIFURCATIONS
    if Dyn_Data.solution_types{1,solution_num}.orbit_type == "free"
        if ~PLOT_BB_SN
            BIFURCATION_SIZE(4) = 0;
        end
    end

    bifurcations = Dyn_Data.bifurcations{1,solution_num};
    bifurcation_types = fields(bifurcations);
    num_bifurcation_types = size(bifurcation_types,1);

    stability = Dyn_Data.stability{1,solution_num};
    [index_ranges,range_stability] = get_stability_boundaries(stability<STABILITY_LIMIT,bifurcations);
    num_sections = size(range_stability,2);

    bifurcation_plot_settings = cell(num_bifurcation_types,1);
    for iType = 1:num_bifurcation_types
        type_plot_settings = {"LineStyle","none","LineWidth",LINE_WIDTH,"Color",line_colour,...
            "Marker",BIFURCATION_MARKER(iType),"MarkerSize",BIFURCATION_SIZE(iType)};
        switch iType
            case 1
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = "w";
                type_plot_settings{end+1} = "MarkerFaceColor";
                type_plot_settings{end+1} = line_colour;
            case 2
                type_plot_settings{end+1} = "MarkerFaceColor";
                type_plot_settings{end+1} = "w";
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = line_colour; %#ok<*AGROW>
            case 3
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = "w";
                type_plot_settings{end+1} = "MarkerFaceColor";
                type_plot_settings{end+1} = line_colour;
        end
        bifurcation_plot_settings{iType,1} =  type_plot_settings;
    end
end

if PLOT_SPECIAL_POINT
    point_index = Dyn_Data.get_special_point(solution_num,"X");
    special_point_plot_settings = bifurcation_plot_settings{1,1}; 
    num_settings = size(special_point_plot_settings,2);
    for iSetting = 1:num_settings
        setting = special_point_plot_settings{1,iSetting};
        if ~isstring(setting)
            continue
        end
        if any(setting == ["Color","MarkerFaceColor"])
            special_point_plot_settings{1,iSetting+1} = get_plot_colours(3);
        end
        if setting == "MarkerSize"
            special_point_plot_settings{1,iSetting+1} = special_point_plot_settings{1,iSetting+1};
        end
        if setting == "LineWidth"
            special_point_plot_settings{1,iSetting+1} = special_point_plot_settings{1,iSetting+1}-0.5;
        end
    end
end

if isempty(Dyn_Data.periodicity_error)
    PLOT_PERIODICITY = 0;
else
    periodicity_error = Dyn_Data.periodicity_error{1,solution_num};
    if isempty(periodicity_error)
        PLOT_PERIODICITY = 0;
    end
end



switch type
    case {"energy","physical amplitude"}
        switch type
            case "energy"
                energy = Dyn_Data.energy{1,solution_num};
            case "physical amplitude"
                energy = Dyn_Data.additional_dynamic_output{1,solution_num};
        end
        

        if isempty(ax)
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

        if PLOT_SPECIAL_POINT
            p = plot(ax,frequency(point_index),energy(point_index),special_point_plot_settings{:});

            data_tip_row_id = dataTipTextRow("ID",orbit_ids(point_index));
            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

            data_tip_row_stab = dataTipTextRow("Stability",stability(point_index));
            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row_stab;

            if PLOT_PERIODICITY
                data_tip_row = dataTipTextRow("Periodicity",periodicity_error(point_index));
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
        Solution_Type = Dyn_Data.solution_types{1,solution_num};
        num_modes = size(amplitude,1);

        if Solution_Type.model_type == "fom"
            r_modes = 1:num_modes;
        else
            r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
        end


        if isempty(ax)
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