function ax = compare_validation(Dyn_Data,type,solution_num,h_modes,varargin)
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%



if isempty(ax)
    fig = figure;
    ax = axes(fig);
end

if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end

if isstring(h_modes) && h_modes == "all"
    h_modes = 1:Dyn_Data.Dynamic_Model.Model.Static_Options.num_validation_modes; 
end

if isstring(solution_num) && solution_num == "last"
    solution_num = Dyn_Data.num_solutions;
end

r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
r_modes = load_data(r_modes);
L_modes = setdiff(h_modes,r_modes);

num_L_modes = length(L_modes);
num_r_modes = length(r_modes);
num_colours = size(get_plot_colours,1);
log_axis = 0;

num_sols = size(solution_num,2);
min_freq = inf;
max_freq = 0;
for iSol = 1:num_sols
    Sol = Dyn_Data.load_solution(solution_num(iSol));
    min_freq = min(min_freq,min(Sol.frequency));
    max_freq = max(max_freq,max(Sol.frequency));
end

for iMode = 1:num_L_modes
    % Dyn_Data.Additional_Output.output = "none";

    switch type
        case "energy"
            plot_backbone = iMode == 1;
        case "amplitude"
            plot_backbone = 0;
        case "force amplitude"
            plot_backbone = 1;
        case "mean error"
            plot_backbone = 1;
        case "validation error"
            plot_backbone = 0;
            log_axis = 1;
    end

    for iSol = 1:num_sols
        [Dyn_Data,Validated_BB_Sol] = Dyn_Data.validate_solution(solution_num(iSol),L_modes(iMode));
        colour_num = mod(iMode-1,num_colours)+1;


        L_mode_index = L_modes(iMode) - (nnz(r_modes < L_modes(iMode)));
        mode_frequency = sqrt(Validated_BB_Sol.low_frequency_eigenvalues(L_mode_index));
        mode_details = sprintf("Last mode: %u - %.2g rad/s - [%.1fx - %.1fx]",[L_modes(iMode),mode_frequency,mode_frequency/max_freq,mode_frequency/min_freq]);

        ax = plot_h_predicition(Dyn_Data,type,solution_num(iSol),"axes",ax,"colour",colour_num,"backbone",plot_backbone);
        if log_axis
            ax.YScale = "log";
        end
        title(ax,mode_details,"FontWeight","normal")

        drawnow
    end
end
% plot_backbone(Dyn_Data,type,solution_num,"axes",ax,"colour",0);
end