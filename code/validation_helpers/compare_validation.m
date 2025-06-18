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
    T = tiledlayout("vertical");
    T.TileSpacing = "none";
    T.Padding = "tight";
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

num_sols = size(solution_num,2);
min_freq = inf;
max_freq = 0;
for iSol = 1:num_sols
    Sol = Dyn_Data.load_solution(solution_num(iSol));
    min_freq = min(min_freq,min(Sol.frequency));
    max_freq = max(max_freq,max(Sol.frequency));
end




num_outputs = numel(type);

backbone_plot = zeros(1,num_outputs);
for iOutput = 1:num_outputs
    nexttile
    ax = gca();
    box(ax,"on")
    switch type(iOutput)
        case "energy"
            backbone_plot(iOutput) = 1;
        case "amplitude"
        case "force amplitude"
            backbone_plot(iOutput) = 1;
        case "mean error"
            backbone_plot(iOutput) = 1;
        case "validation error"
            ax.YScale = "log";
            ylim(ax,[1e-5,1])
            ylabel("\epsilon","Interpreter","tex")
            validation_error_labels = ax.YTickLabel;
       
        case "physical amplitude"
            backbone_plot(iOutput) = 1;
    end

    if backbone_plot(iOutput)
        ax = compare_solutions(type(iOutput),Dyn_Data,solution_num,"axes",ax,"legend",0);
        bb_lines = findobj(ax,"Type","line");
        set(bb_lines,"Tag","backbone");
    end



    if iOutput < num_outputs
        ax.XTickLabel = repmat("",size(ax.XTickLabel));
    end
    if iOutput > 1
        ax.YTickLabel{end} = "";
    end
end
ax = flip(findobj(fig,"Type","axes"));


ax(type == "validation error").YLim = [1e-5,1.1];
% drawnow


for iMode = 1:num_L_modes
    % Dyn_Data.Additional_Output.output = "none";



    for iSol = 1:num_sols
        [Dyn_Data,Validated_BB_Sol] = Dyn_Data.validate_solution(solution_num(iSol),L_modes(iMode));
        colour_num = mod(iMode-1,num_colours)+1;


        L_mode_index = L_modes(iMode) - (nnz(r_modes < L_modes(iMode)));
        mode_frequency = sqrt(Validated_BB_Sol.low_frequency_eigenvalues(L_mode_index));
        mode_details = sprintf("Last mode: %u - %.2g rad/s - [%.1fx - %.1fx]",[L_modes(iMode),mode_frequency,mode_frequency/max_freq,mode_frequency/min_freq]);
        
        for iOutput = 1:num_outputs
            plot_h_predicition(Dyn_Data,type(iOutput),solution_num(iSol),"axes",ax(iOutput),"colour",colour_num,"backbone",0);
            if iOutput < num_outputs
                ax(iOutput).XTickLabel = repmat("",size(ax(iOutput).XTickLabel));
            end
            if iOutput > 1
                ax(iOutput).YTickLabelMode = "auto";
                ax(iOutput).YTickLabel{end} = "";
            end
            if backbone_plot(iOutput)
                bb_lines = findobj(ax(iOutput),"Tag","backbone");
                uistack(bb_lines,"top")
            end
            ax(iOutput).YLim = ax(iOutput).YLim;

            lines = findobj(ax(iOutput),"Type","line");
            num_lines = numel(lines);
            for iLine = 1:num_lines
                if min(lines(iLine).YData) < ax(iOutput).YLim(1)
                    ax(iOutput).YLim(1) = min(lines(iLine).YData);
                end
            end
        end
        title(ax(1),mode_details,"FontWeight","normal")
        drawnow
    end
end
% plot_backbone(Dyn_Data,type,solution_num,"axes",ax,"colour",0);
end