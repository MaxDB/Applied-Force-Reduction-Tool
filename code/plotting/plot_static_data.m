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
colour = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case "outputs"
            outputs = keyword_values{arg_counter};
        case "plot seps"
            plot_seps = keyword_values{arg_counter};
        case {"colour","color"}
            colour = keyword_values{arg_counter};
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
is_matrix = 0;
switch type
    case {"displacement","physical_displacement"}
        y_data = Static_Data.get_dataset_values("physical_displacement");
        y_origin = zeros(size(y_data,1),1);
        y_label = "x";
    case {"force","restoring_force"}
        y_data = Static_Data.get_dataset_values("restoring_force");
        y_origin = zeros(size(y_data,1),1);
        y_label = "f";
    case {"potential","potential_energy","energy"}
        y_data = Static_Data.get_dataset_values("potential_energy");
        y_origin = zeros(size(y_data,1),1);
        y_label = "V";
    case{"force-displacement"}
        y_data = Static_Data.get_dataset_values("physical_displacement");
        y_origin = zeros(size(y_data,1),1);
        y_label = "x";
        x_data = Static_Data.get_dataset_values("restoring_force");
        x_label = "f";
        x_origin = repmat(R_ORIGIN,num_r_modes,1); %#ok<RPMT0>
    case {"force-energy"}
        y_data = Static_Data.get_dataset_values("potential_energy");
        y_origin = zeros(size(y_data,1),1);
        y_label = "V";
        x_data = Static_Data.get_dataset_values("restoring_force");
        x_label = "f";
        x_origin = repmat(R_ORIGIN,num_r_modes,1); %#ok<RPMT0>
    case "h_stiffness"
        y_data = Static_Data.get_dataset_values("low_frequency_stiffness");
        is_matrix = 1;
        mat_size = size(y_data);
        vec_size = {mat_size(1)*mat_size(2),mat_size(3)};
        y_data = reshape(y_data,vec_size{:});

        Model = Static_Data.Model;
        [~,~,v_modeshapes] = Static_Data.get_current_h_data;
        h_transform = v_modeshapes'*Model.mass;

        H_0 = h_transform*(Model.stiffness\(h_transform'));
        stiffness_0 = eye(size(v_modeshapes,2))/H_0;

        y_origin = reshape(stiffness_0,vec_size{1},1);
        y_label = "D";
    case "h_displacement_gradient"
        y_data = Static_Data.get_dataset_values("low_frequency_coupling_gradient");
        is_matrix = 1;
        mat_size = size(y_data);
        vec_size = {mat_size(1)*mat_size(2),mat_size(3)};
        y_data = reshape(y_data,vec_size{:});

        Model = Static_Data.Model;
        [~,~,v_modeshapes] = Static_Data.get_current_h_data;
        h_transform = v_modeshapes'*Model.mass;
        X_0 = (Model.stiffness\(h_transform'));
        H_0 = h_transform*X_0;
        G_0 = X_0/H_0;

        y_origin = reshape(G_0,vec_size{1},1);
        y_label = "G";
    case "perturbation"
        y_data = Static_Data.get_dataset_values("perturbation_displacement");
        is_matrix = 1;
        mat_size = size(y_data);
        vec_size = {mat_size(1)*mat_size(2),mat_size(3)};
        y_data = reshape(y_data,vec_size{:});
        lambda = Static_Data.perturbation_scale_factor;
        Model = Static_Data.Model;
        
        v_modeshapes = [Model.reduced_eigenvectors,Model.low_frequency_eigenvectors];
        x_0 = Model.stiffness\(Model.mass*v_modeshapes.*lambda);
        y_origin = reshape(x_0,numel(x_0),1);
        y_label = "\tilde{\mathbf x}^*";
    otherwise
        error("Plotting for '" + type + "' is not supported")
end


if isempty(outputs)
    total_outputs = size(y_data,1);
    if total_outputs > MAX_OUTPUTS
        outputs = zeros(1,MAX_OUTPUTS);
        while length(unique(outputs)) < MAX_OUTPUTS
            outputs = randi(total_outputs,1,MAX_OUTPUTS);
        end
        outputs = unique(outputs);
    else
        outputs = 1:total_outputs;
    end
elseif is_matrix
    % outputs = sub2ind([mat_size(2),mat_size(1)],outputs(:,2),outputs(:,1));
    outputs = sub2ind([mat_size(1),mat_size(2)],outputs(:,1),outputs(:,2));
end

    



y_data = y_data(outputs,:);
y_origin = y_origin(outputs);
num_outputs = length(outputs);

%-------------------------------------------------------------------------%
%data ordering
sep_id = Static_Data.static_equilibrium_path_id;
num_seps = max(sep_id);

if plot_seps
    marker_type = ".-";
    sep_ordering = Static_Data.get_dataset_values("restoring_force");
else
    marker_type = ".";
    sep_ordering = ones(1,num_points);
end

%-------------------------------------------------------------------------%
%plot style

plot_settings = {marker_type,"MarkerSize",MARKER_SIZE,"LineWidth",LINE_WIDTH};
if ~isempty(colour)
    plot_settings = [plot_settings,{"Color",colour}];
end
origin_plot_settings = [plot_settings,{"MarkerEdgeColor",[0,0,0]}];

%-------------------------------------------------------------------------%
%plotting 
if isa(ax,"matlab.graphics.axis.Axes")
    ax = {ax};
    
end

axes_exist = ~isempty(ax);
if ~axes_exist
    figure
    tiledlayout("flow")
    ax = cell(num_outputs,1);
end

if num_r_modes > 2
    warning("Would require a " + (num_r_modes + 1) + "D plot")
    return
end

for iOutput = 1:num_outputs
    if ~axes_exist
        s1 = nexttile;
    else
        s1 = ax{iOutput};
    end
    box(s1,"on")
    hold(s1,"on")

    for iSep = 1:num_seps
        sep_span = sep_id == iSep;
        x_sep = x_data(:,sep_span);
        y_sep = y_data(:,sep_span);
        force_sep = sep_ordering(:,sep_span);
        [~,sep_order] = sort(abs(force_sep(1,:)),"ascend");

        x_plot = [x_origin,x_sep(:,sep_order)];
        y_plot = [y_origin(iOutput),y_sep(iOutput,sep_order)];



        switch num_r_modes
            case 1
                plot(s1,x_plot,y_plot,plot_settings{:})
                plot(s1,x_origin,y_origin(iOutput),origin_plot_settings{:})

                xlabel(s1,"$" + x_label + "_{" + r_modes + "}$","Interpreter","latex")
                if ~is_matrix
                    y_output_label = y_label + "_{" + outputs(iOutput) + "}";
                else
                    [output_row,output_col] = ind2sub([mat_size(1),mat_size(2)],outputs(iOutput));
                    y_output_label = y_label + "_{(" + output_row + "," + output_col + ")}";
                end
                ylabel(s1,"$" + y_output_label + "$","Interpreter","latex")
            case 2
                plot3(s1,x_plot(1,:),x_plot(2,:),y_plot,plot_settings{:})
                plot3(s1,x_origin(1,:),x_origin(2,:),y_origin(iOutput),origin_plot_settings{:})

                xlabel(s1,"$" + x_label + "_{" + r_modes(1) + "}$","Interpreter","latex")
                ylabel(s1,"$" + x_label + "_{" + r_modes(2) + "}$","Interpreter","latex")
                if ~is_matrix
                    z_output_label = y_label + "_{" + outputs(iOutput) + "}";
                else
                    output_row = ceil(outputs(iOutput)/mat_size(2));
                    output_col = outputs(iOutput) - (output_row-1)*mat_size(2);
                    z_output_label = y_label + "_{(" + output_row + "," + output_col + ")}";
                end
                zlabel(s1,"$" + z_output_label + "$","Interpreter","latex")
        end
    end

    hold(s1,"off")
    ax{iOutput} = s1; %#ok<AGROW>
end



end