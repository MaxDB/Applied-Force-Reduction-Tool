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
validated = 1;

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
        case {"validation"}
            validated = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
line_colour = get_plot_colours(colour_num);
line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",line_colour,"DisplayName","r"};

if validated
    validated_line_colour = get_plot_colours(colour_num+1);
    validated_line_plot_settings = {"LineWidth",LINE_WIDTH,"Color",validated_line_colour,"DisplayName","h"};
end


if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end
[Orbit,Validated_Orbit] = Dyn_Data.get_orbit(solution_num,orbit_num,validated);
Rom = Dyn_Data.Dynamic_Model;
Model = Rom.Model;

Sol = Dyn_Data.load_solution(solution_num);
Sol_Type = Dyn_Data.solution_types{1,solution_num};

known_modes = Model.reduced_modes;
known_eval = Model.reduced_eigenvalues;
known_evec = Model.reduced_eigenvectors.load();
if validated && Sol_Type.validated
    Sol_v = Dyn_Data.load_solution(solution_num,"validation");

    lf_modes = Model.low_frequency_modes;
    known_modes = [known_modes,lf_modes];
   
    orbit_lf_modes = Sol_v.validation_modes;
    r_modes = Model.reduced_modes;
    orbit_h_modes = [r_modes,orbit_lf_modes];

    
    lf_eval = Model.low_frequency_eigenvalues;
    known_eval = [known_eval;lf_eval];
    
    lf_evec = Model.low_frequency_eigenvectors.load();
    known_evec = [known_evec,lf_evec];

  
end
modal_transform = known_evec'*Model.mass;

time = Orbit.tbp';
state = Orbit.xbp';
state_size = size(state,1);

num_modes = state_size/2;
disp_index = 1:num_modes;
vel_index = disp_index + num_modes;


if normalise
    period = Orbit.T;
    time = time/period;
end

plot_dimension = length(type);
num_time_points = size(time,2);

orbit_plot = zeros(plot_dimension,num_time_points);
labels = strings(plot_dimension,1);
for iType = 1:plot_dimension
    data_id = split(type(iType),"-");
    plot_type = data_id(1); %physical vs modal vs time etc.
    plot_part = data_id(2); % displacement vs velocity
    data_index = double(data_id(3)); %specific quantity: mode 2 etc..

    label = "";
    switch plot_type
        case "r"
            plotted_state = state;
            plotted_output_size = num_modes;
            label = label + plot_type;
        case "h"
            plotted_state = [Validated_Orbit.h;Validated_Orbit.h_dot];
            plotted_output_size = size(Validated_Orbit.h,1);
            label = label + plot_type;
        case "q"
            x = Rom.expand(state(disp_index,:));
            x_vel = Rom.expand_velocity(state(disp_index,:),state(vel_index,:));
            
            mode_map = known_modes == data_index;
            q_transform = modal_transform(mode_map,:);
            q = q_transform*x;
            q_vel = q_transform*x_vel;

            plotted_state = [q;q_vel];
            plotted_output_size = size(q,1);
            label = label + plot_type;
        case "v"
            x = Rom.expand(state(disp_index,:));
            x_vel = Rom.expand_velocity(state(disp_index,:),state(vel_index,:));
            
            mode_map = ismember(known_modes,orbit_h_modes);
            v_transform = modal_transform(mode_map,:);
            q = v_transform*x;
            q_vel = v_transform*x_vel;

            plotted_state = [q;q_vel] + [Validated_Orbit.h;Validated_Orbit.h_dot];
            plotted_output_size = size(q,1);
            label = label + plot_type;
        case "x"
            %plotted_state = expands...
        case "t"
            data_id(2:3) = ["0";"1"];
            plotted_value_type = time;
            label = label + plot_type;
    end

    switch plot_part
        case "d"
            plotted_value_type = plotted_state(1:plotted_output_size,:);
        case "v"
            plotted_value_type = plotted_state((plotted_output_size + 1):(2*plotted_output_size),:);
            label = "\dot{" + label + "}";
    end

    
    switch plot_type
        case {"r","x","h","v","q"}
            label = label + "_{" + data_index + "}";
    end

    switch plot_type
        case {"q"}
            plotted_value = plotted_value_type(1,:);
        otherwise
            plotted_value = plotted_value_type(data_index,:);
    end
    
    orbit_plot(iType,:) = plotted_value; 
    labels(iType) = "$" + label + "$";
end

%-------------------
if isempty(ax)
    fig = figure;
    ax = axes(fig);
    box(ax,"on")
end
hold(ax,"on")
switch plot_dimension
    case 2
        p = plot(ax,orbit_plot(1,:),orbit_plot(2,:),line_plot_settings{:});
        xlabel(ax,labels(1),"Interpreter","latex")
        ylabel(ax,labels(2),"Interpreter","latex")
    case 3
        p = plot3(ax,orbit_plot(1,:),orbit_plot(2,:),orbit_plot(3,:),line_plot_settings{:});
        xlabel(ax,labels(1),"Interpreter","latex")
        ylabel(ax,labels(2),"Interpreter","latex")
        zlabel(ax,labels(3),"Interpreter","latex")
end
hold(ax,"off")

%add cursor data
id_data = "[" + solution_num + "," + orbit_num + "]";
data_tip_row = dataTipTextRow("ID",repmat(id_data,1,num_time_points));
p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;

if normalise 
    time_label = "t/T";
else
    time_label = "Time";
end
data_tip_row = dataTipTextRow(time_label,time);
p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;


data_tip_row = dataTipTextRow("Period",repmat(Orbit.T,1,num_time_points));
p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;


r_modes = Model.reduced_modes;
reduction_basis = "\{" + join(string(r_modes),",") + "\}";
basis_label = reduction_basis;
if plot_type == "v"
    basis_label = basis_label + " : " + "\{" + join(string(orbit_h_modes),",") + "\}";
end
data_tip_row = dataTipTextRow("Basis",repmat(basis_label,1,num_time_points));
p.DataTipTemplate.DataTipRows(end+1) = data_tip_row;
% switch type
%     case {"displacement","phase"}
%         displacement = orbit.xbp';
%         num_r_modes = size(displacement,1)/2;
% 
%         velocity = displacement((num_r_modes+1):(2*num_r_modes),:);
%         displacement = displacement(1:num_r_modes,:);
% 
% 
%         Solution_Type = Dyn_Data.solution_types{1,solution_num};
% 
% 
% 
%         if Solution_Type.model_type == "fom"
%             r_modes = 1:num_r_modes;
%         else
%             r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
%         end
% 
%         if validated
%             Validated_Solution = Dyn_Data.load_solution(solution_num,"validation");
%             L_modes = Validated_Solution.validation_modes;
% 
%             physical_displacement = Dyn_Data.Dynamic_Model.expand(displacement);
%             Model = Dyn_Data.Dynamic_Model.Model;
%             mass = Model.mass;
%             % stiffness = Model.stiffness;
%             L_evecs = Dyn_Data.Dynamic_Model.Model.low_frequency_eigenvectors;
% 
%             all_L_modes = 1:max(Dyn_Data.Dynamic_Model.Model.low_frequency_modes);
%             all_L_modes(ismember(all_L_modes,r_modes)) = [];
%             L_map = ismember(all_L_modes,L_modes);
% 
%             L_evecs = L_evecs(:,L_map);
% 
% 
%             % r_evecs = Model.reduced_eigenvectors;
% 
%             g = L_evecs'*mass*physical_displacement;
% 
%             displacement = [displacement;g];
% 
%             h_displacement = validated_orbit.h;
%             h_displacement = displacement + h_displacement;
% 
% 
%             h_modes = [r_modes,L_modes];
%             plot_modes = h_modes;
%         else
%             plot_modes = r_modes;
%         end
% 
%         num_modes = size(plot_modes,2);
% 
%         if exist("h_displacement","var")
%             displacement_all = [displacement;h_displacement];
%         else
%             displacement_all = displacement;
%         end
%         [time_shifted,displacement_shifted] = shift_orbit(time,displacement_all);
%         time = time_shifted;
%         displacement = displacement_shifted(1:num_modes,:);
%         if exist("h_displacement","var")
%             h_displacement = displacement_shifted((num_modes+1):(2*num_modes),:);
%         end
% 
%         switch type
%             case "displacement"
% 
%                 if isempty(ax)
%                     figure
%                     tiledlayout("flow")
%                     for iMode = 1:num_modes
%                         ax{iMode,1} = nexttile; %#ok<*AGROW>
%                         ax{iMode,2} = plot_modes(iMode);
%                     end
%                     plotted_modes = plot_modes;
%                 else
%                     num_axes = size(ax,1);
%                     %check if all required r_modes are present]
%                     for iAx = 1:num_axes
%                         plotted_modes(iAx,1) = ax{iAx,2};
%                     end
% 
%                     neglected_modes = setdiff(r_modes,plotted_modes);
%                     num_neglected_modes = length(neglected_modes);
%                     for iMode = 1:num_neglected_modes
%                         ax{num_axes+iMode,1} = nexttile;
%                         ax{num_axes+iMode,2} = neglected_modes(iMode);
%                         plotted_modes(num_axes+iMode) = neglected_modes(iMode);
%                     end
%                 end
% 
% 
% 
%                 for iMode = 1:num_modes
%                     mode = plot_modes(iMode);
%                     ax_id = plotted_modes == mode;
%                     box(ax{ax_id,1},"on")
% 
% 
%                     hold(ax{ax_id,1},"on")
%                     plot(ax{ax_id,1},time,displacement(iMode,:),line_plot_settings{:});
%                     if validated
%                         plot(ax{ax_id,1},time,h_displacement(iMode,:),validated_line_plot_settings{:});
%                     end
%                     hold(ax{ax_id,1},"off")
% 
%                     xlabel(ax{ax_id,1},"Time (s)")
%                     ylabel(ax{ax_id,1},"r_{" + plot_modes(iMode) + "}")
%                     xlim(ax{ax_id},[min(time),max(time)])
%                     legend(ax{ax_id})
%                 end
%                 title(ax{1,1},"(" + solution_num + "," + orbit_num + ")")
%             case "phase"
%                 if isempty(ax)
%                     figure;
%                     ax = gca;
%                     box(ax,"on")
%                     title(ax,"(" + solution_num + "," + orbit_num + ")")
% 
%                     if validated
%                         xlabel("r") 
%                         ylabel("s") 
%                     else
%                         xlabel("r") 
%                         ylabel("r_{dot}") 
%                     end
%                 end
% 
%                 if num_r_modes == 1
%                     displacement = [displacement;velocity];
%                 end
% 
%                 if validated
%                     hold(ax,"on")
%                     plot(ax,displacement(1,:),displacement(2,:),line_plot_settings{:})
%                     plot(ax,h_displacement(1,:),h_displacement(2,:),validated_line_plot_settings{:})
%                     hold(ax,"off")
%                     xlabel("r")
%                     ylabel("s")
%                 else
%                     plot(ax,displacement(1,:),displacement(2,:),line_plot_settings{:})
%                 end
% 
%         end
%     case "physical displacement"
%         displacement = orbit.xbp';
%         num_modes = size(displacement,1)/2;
%         displacement = displacement(1:num_modes,:);
% 
%         Rom = Dyn_Data.Dynamic_Model;
%         node_map = Rom.Model.node_mapping;
%         dof_bc = Dyn_Data.Additional_Output.dof;
%         dof = node_map(node_map(:,1) == dof_bc,2);
% 
%         x_phy = Rom.expand(displacement);
%         x_dof = x_phy(dof,:);
% 
%         if isempty(ax)
%             figure;
%             ax = gca;
%             box(ax,"on")
%             ylabel(ax,"x_{" + dof_bc + "} (m)")
%             if normalise
%                 xlabel(ax,"Normalised Time")
%                 xlim(ax,[0,1])
%             else
%                 xlabel(ax,"Time (s)")
%             end
%         end
% 
%         if orbit_shift
%             [time,x_dof] = shift_orbit(time,x_dof);
%         end
% 
%         hold(ax,"on")
%         plot(time,x_dof,line_plot_settings{:})
%         hold(ax,"off")
%     case "max physical displacement"
%         displacement = orbit.xbp';
%         num_modes = size(displacement,1)/2;
%         displacement = displacement(1:num_modes,:);
% 
%         Rom = Dyn_Data.Dynamic_Model;
%         node_map = Rom.Model.node_mapping;
%         dof_bc = Dyn_Data.Additional_Output.dof;
%         dof = node_map(node_map(:,1) == dof_bc,2);
% 
%         x_phy = Rom.expand(displacement);
%         x_dof = x_phy(dof,:);
%         [~,max_disp_index] = max(abs(x_dof));
%         max_disp = x_phy(:,max_disp_index);
% 
%         Model = Dyn_Data.Dynamic_Model.Model;
%         plot_fe_mesh(Model,max_disp)
% end

end