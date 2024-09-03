function compare_validation(Dyn_Data,type,solution_num,L_modes)
ax = [];

num_L_modes = length(L_modes);
num_colours = size(get_plot_colours,1);

for iMode = 1:num_L_modes
    Dyn_Data = Dyn_Data.validate_solution(solution_num,L_modes(iMode));
    colour_num = mod(iMode-1,num_colours)+1;

    switch type
        case "energy"
            plot_backbone = iMode == num_L_modes;
        case "amplitude"
            plot_backbone = 1;
        case "force amplitude"
            plot_backbone = 1;
    end
           
    ax = plot_h_predicition(Dyn_Data,type,solution_num,"axes",ax,"colour",colour_num,"backbone",plot_backbone);
end
% plot_backbone(Dyn_Data,type,solution_num,"axes",ax,"colour",0);
end