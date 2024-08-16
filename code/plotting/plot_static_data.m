function ax = plot_static_data(type,Static_Data,varargin)
MAX_OUTPUTS = 4;
R_ORIGIN = 0;

LINE_WIDTH = 1;
MARKER_SIZE = 6;

%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
outputs = [];
plot_seps = 1;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case "outputs"
            outputs = keyword_values{arg_counter};
        case {"plot seps"}
            plot_seps = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
if isstring(Static_Data)
    Static_Data = load_static_data(Static_Data);
end

%-------------------------------------------------------------------------%
%input data

x_data = Static_Data.get_dataset_values("reduced_displacement");
x_label = "r";

r_modes = Static_Data.Model.reduced_modes;
num_r_modes = size(x_data,1);
num_points = size(x_data,2);

x_origin = repmat(R_ORIGIN,num_r_modes,1); %#ok<RPMT0>
%-------------------------------------------------------------------------%
%output data
switch type
    case {"displacement","physical_displacement"}
        y_data = Static_Data.get_dataset_values("physical_displacement");
        y_origin = 0;
        y_label = "x";
    case {"force","restoring_force"}
        y_data = Static_Data.get_dataset_values("restoring_force");
        y_origin = 0;
        y_label = "f";
    case {"potential","potential_energy","energy"}
        y_data = Static_Data.get_dataset_values("potential_energy");
        y_origin = 0;
        y_label = "V";
    case{"force-displacement"}
        y_data = Static_Data.get_dataset_values("physical_displacement");
        y_origin = 0;
        y_label = "x";

        x_data = Static_Data.get_dataset_values("restoring_force");
        x_label = "f";

        x_origin = repmat(R_ORIGIN,num_r_modes,1); %#ok<RPMT0>
    otherwise
        error("Plotting for '" + type + "' is not supported")
end


switch ndims(y_data)
    case {1,2}
        if isempty(outputs)
            total_outputs = size(y_data,1);
            if total_outputs > MAX_OUTPUTS
                outputs = randi(total_outputs,1,MAX_OUTPUTS);
            else
                outputs = 1:total_outputs;
            end
        end

        y_data = y_data(outputs,:);
    case {3}
        y_data = y_data(outputs,:,:);
end
num_outputs = length(outputs);

%-------------------------------------------------------------------------%
%data ordering
sep_id = Static_Data.static_equilibrium_path_id;
num_seps = max(sep_id);

if plot_seps
    marker_type = ".-";
    sep_ordering = Static_Data.get_dataset_values("restoring_force");
else
    marker_type = "x";
    sep_ordering = ones(1,num_points);
end

%-------------------------------------------------------------------------%
%plot style

plot_settings = {marker_type,"MarkerSize",MARKER_SIZE,"LineWidth",LINE_WIDTH};
origin_plot_settings = [plot_settings,{"MarkerEdgeColor",[0,0,0]}];

%-------------------------------------------------------------------------%
%plotting 
if isa(ax,"matlab.graphics.axis.Axes")
    ax = {ax};
end

if isempty(ax)
    figure
    tiledlayout("flow")
    ax = cell(num_outputs,1);
end

if num_r_modes > 2
    warning("Would require a " + (num_r_modes + 1) + "D plot")
    return
end

for iOutput = 1:num_outputs
    s1 = nexttile;
    box on
    hold(s1,"on")

    for iSep = 1:num_seps
        sep_span = sep_id == iSep;
        x_sep = x_data(:,sep_span);
        y_sep = y_data(:,sep_span);
        force_sep = sep_ordering(:,sep_span);
        [~,sep_order] = sort(abs(force_sep(1,:)),"ascend");

        x_plot = [x_origin,x_sep(:,sep_order)];
        y_plot = [y_origin,y_sep(iOutput,sep_order)];



        switch num_r_modes
            case 1
                plot(x_plot,y_plot,plot_settings{:})
                plot(x_origin,y_origin,origin_plot_settings{:})

                xlabel(x_label + "_{" + r_modes + "}")
                if num_outputs == 1
                    ylabel(y_label)
                else
                    ylabel(y_label + "_{" + outputs(iOutput) + "}")
                end
            case 2
                plot3(x_plot(1,:),x_plot(2,:),y_plot,plot_settings{:})
                plot3(x_origin(1,:),x_origin(2,:),y_origin,origin_plot_settings{:})

                xlabel(x_label + "_{" + r_modes(1) + "}")
                ylabel(x_label + "_{" + r_modes(2) + "}")
                zlabel(y_label + "_{" + outputs(iOutput) + "}")

        end
    end

    hold(s1,"off")
    ax{iOutput} = s1; %#ok<AGROW>
end



end