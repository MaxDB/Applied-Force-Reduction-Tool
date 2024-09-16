function ax = plot_backbone(Dyn_Data,type,solution_num,varargin)
PLOT_BIFURCATIONS = 1;
PLOT_SPECIAL_POINT = 0;

STABILITY_LIMIT = 1.005;

LINE_STYLE = [":","-"]; %[unstable,stable]
LINE_WIDTH = 1.5;
BIFURCATION_MARKER = ["o","o","^","x"];
BIFURCATION_TYPE = ["BP","PD","NS","SN"];

BB_BIFURCATION_SIZE = [4,4,4,0];
FRF_BIFURCATION_SIZE = [4,4,4,5];
SPECIAL_POINT_SIZE = 6;

FE_DATA_MARKER = "*";
FE_DATA_COLOUR = [0,0,0; 1,0,0];
FE_DATA_MARKER_SIZE = 6;
FE_PERIODICITY_LIMIT = 0.01;
SHOW_UNCONVERGED = 1;

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
Solution = Dyn_Data.load_solution(solution_num);
frequency = Solution.frequency;
num_orbits = Solution.num_orbits;
orbit_labels = 1:num_orbits;
orbit_ids = "(" + solution_num + "," + orbit_labels + ")";

line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",line_colour};

extra_data = [];
if Solution.Solution_Type.orbit_type == "forced"
    if isfield(Solution.Solution_Type,"amplitude")
        extra_data = Solution.Solution_Type.amplitude;
        extra_data_name = "Amplitude";
    end
    
end

if PLOT_BIFURCATIONS
    switch Solution.Solution_Type.orbit_type
        case "free"
            bifurcation_size = BB_BIFURCATION_SIZE;
        case "forced"
            bifurcation_size = FRF_BIFURCATION_SIZE;
    end

    bifurcations = Solution.bifurcations;
    bifurcation_types = fields(bifurcations);
    num_bifurcation_types = size(bifurcation_types,1);

    stability = Solution.stability;
    [index_ranges,range_stability] = get_stability_boundaries(stability<STABILITY_LIMIT,bifurcations);
    num_sections = size(range_stability,2);

    bifurcation_plot_settings = cell(num_bifurcation_types,1);
    for iType = 1:num_bifurcation_types
        type_plot_settings = {"LineStyle","none","LineWidth",LINE_WIDTH,"Color",line_colour,...
            "Marker",BIFURCATION_MARKER(iType),"MarkerSize",bifurcation_size(iType)};
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
                type_plot_settings{end+1} = "MarkerFaceColor";
                type_plot_settings{end+1} = "w";
                type_plot_settings{end+1} = "MarkerEdgeColor";
                type_plot_settings{end+1} = line_colour;
        end
        bifurcation_plot_settings{iType,1} =  type_plot_settings;
    end
end

if PLOT_SPECIAL_POINT == 1
    plot_special_point = ~isempty(Dyn_Data.Additional_Output);
else
    plot_special_point = PLOT_SPECIAL_POINT;
end

if plot_special_point
    point_index = Dyn_Data.get_special_point(solution_num,"X");
    if isempty(point_index)
        plot_special_point = 0;
    end
    special_point_plot_settings = bifurcation_plot_settings{1,1}; 
    marker_size_index = find(cellfun(@(iSetting) isequal(iSetting,"MarkerSize"),special_point_plot_settings));
    special_point_plot_settings{1,marker_size_index+1} = SPECIAL_POINT_SIZE;
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

Periodicity_Ouput = Dyn_Data.load_solution(solution_num,"periodicity");
plot_periodicity = class(Periodicity_Ouput) == "FE_Orbit_Output";
if plot_periodicity
    periodicity_error = nan(1,num_orbits);
    periodicity_error(Periodicity_Ouput.orbit_labels) = Periodicity_Ouput.fe_output;
end

FE_Output = Dyn_Data.load_solution(solution_num,"forced_response");
plot_fe_output = class(FE_Output) == "FE_Orbit_Output";
if plot_fe_output
    FE_Data = FE_Output.fe_output;
    converged_fe = FE_Data.periodicity < FE_PERIODICITY_LIMIT;
    fe_data_plot_settings = {"Marker",FE_DATA_MARKER,"MarkerSize",FE_DATA_MARKER_SIZE,"LineStyle","none"};
end


switch type
    case {"energy","physical amplitude","stress"}
        switch type
            case "energy"
                energy = Solution.energy;
            case "physical amplitude"
                energy = Solution.additional_dynamic_output;
            case "stress"
                % stress_dof = Dyn_Data.Additional_Output.dof;
                % stress_node = 1 + (stress_dof - mod(stress_dof,6))/6;
                energy = max(Solution.max_displacement_stress,[],1);
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

                if plot_periodicity
                    data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end

                if ~isempty(extra_data)
                    data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(index_range)));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end

            end

            for iType = 1:num_bifurcation_types
                bifurcation_type = bifurcation_types{iType,1};
                bifurcation_index = bifurcations.(bifurcation_type);
                if isempty(bifurcation_index)
                    continue
                end
                if bifurcation_size(iType) == 0
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
                

                if plot_periodicity
                    data_tip_row = dataTipTextRow("Periodicity",periodicity_error(bifurcation_index));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end

                if ~isempty(extra_data)
                    data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(bifurcation_index)));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end
            end

        else
            p = plot(ax,frequency,energy,'LineStyle',LINE_STYLE(2),line_plot_settings{:});
            data_tip_row = dataTipTextRow("ID",orbit_ids);
            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

            if plot_periodicity
                data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
            end

            if ~isempty(extra_data)
                data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(index_range)));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
            end
        end

        if plot_special_point
            p = plot(ax,frequency(point_index),energy(point_index),special_point_plot_settings{:});

            data_tip_row_id = dataTipTextRow("ID",orbit_ids(point_index));
            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row_id;

            data_tip_row_stab = dataTipTextRow("Stability",stability(point_index));
            p.DataTipTemplate.DataTipRows(end+1) = data_tip_row_stab;

            if plot_periodicity
                data_tip_row = dataTipTextRow("Periodicity",periodicity_error(point_index));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
            end

            if ~isempty(extra_data)
                data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(point_index)));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
            end
        end

        if plot_fe_output
            switch type
                case "energy"
                    fe_data_plot = FE_Data.energy;
                case "physical amplitude"
                    fe_data_plot = FE_Data.additional_dynamic_output;
            end
            fe_plot_index = converged_fe;
            num_fe_plots = SHOW_UNCONVERGED + 1;
            for iFE_plot = 1:num_fe_plots
                if nnz(fe_plot_index) == 0
                    fe_plot_index = ~fe_plot_index;
                    continue
                end
                p = plot(FE_Data.frequency(fe_plot_index ),fe_data_plot(fe_plot_index ),"Color",FE_DATA_COLOUR(iFE_plot,:),fe_data_plot_settings{:});

                data_tip_row = dataTipTextRow("ID",orbit_ids(arrayfun(@(label) find(orbit_labels == label),FE_Output.orbit_labels(fe_plot_index ))));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                
                data_tip_row = dataTipTextRow("Periodicity",FE_Data.periodicity(fe_plot_index));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                data_tip_row = dataTipTextRow("Simulated periods",FE_Data.simulated_periods(fe_plot_index));
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                if ~isempty(extra_data)
                    data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(fe_plot_index)));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end
                fe_plot_index = ~fe_plot_index;
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
       
        amplitude = Solution.amplitude;
        Solution_Type = Solution.Solution_Type;
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

                    if plot_periodicity
                        data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end


                    if ~isempty(extra_data)
                        data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(index_range)));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end


                end

                for iType = 1:num_bifurcation_types
                    bifurcation_type = bifurcation_types{iType,1};
                    bifurcation_index = bifurcations.(bifurcation_type);
                    if isempty(bifurcation_index)
                        continue
                    end
                    if bifurcation_size(iType) == 0
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

                    if plot_periodicity
                        data_tip_row = dataTipTextRow("Periodicity",periodicity_error(bifurcation_index));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end

                    if ~isempty(extra_data)
                        data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(bifurcation_index)));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end
                end
            else
                p = plot(ax{ax_id,1},frequency,amplitude(iMode,:),'LineStyle',LINE_STYLE(2),line_plot_settings{:});
                data_tip_row = dataTipTextRow("ID",orbit_ids);
                p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                if plot_periodicity
                    data_tip_row = dataTipTextRow("Periodicity",periodicity_error(index_range));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                end
            end
            if plot_fe_output

                fe_plot_index = converged_fe;
                num_fe_plots = SHOW_UNCONVERGED + 1;
                for iFE_plot = 1:num_fe_plots
                    if nnz(fe_plot_index) == 0
                        fe_plot_index = ~fe_plot_index;
                        continue
                    end
                    p = plot(ax{ax_id,1},FE_Data.frequency(fe_plot_index ),FE_Data.amplitude(iMode,fe_plot_index),"Color",FE_DATA_COLOUR(iFE_plot,:),fe_data_plot_settings{:});
                    
                    data_tip_row = dataTipTextRow("ID",orbit_ids(arrayfun(@(label) find(orbit_labels == label),FE_Output.orbit_labels(fe_plot_index ))));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    data_tip_row = dataTipTextRow("Periodicity",FE_Data.periodicity(fe_plot_index));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    data_tip_row = dataTipTextRow("Simulated periods",FE_Data.simulated_periods(fe_plot_index));
                    p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

                    if ~isempty(extra_data)
                        data_tip_row = dataTipTextRow(extra_data_name,extra_data*ones(size(fe_plot_index)));
                        p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
                    end

                    fe_plot_index = ~fe_plot_index;
                end
            end


            hold(ax{ax_id,1},"off")

            xlabel(ax{ax_id,1},"Frequency (rad/s)")
            ylabel(ax{ax_id,1},"R_{" + r_modes(iMode) + "}")


        end
end
end