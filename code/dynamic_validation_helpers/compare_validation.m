function compare_validation(type,Dyn_Data,solution_num,L_modes)
switch type
    case "energy"
        figure
        ax = gca;

        num_L_modes = length(L_modes);
        num_colours = size(get_plot_colours,1);
        for iMode = 1:num_L_modes
            Dyn_Data = Dyn_Data.validate_solution(solution_num,L_modes(iMode));
            colour_num = mod(iMode-1,num_colours)+1;
            ax = plot_h_predicition(Dyn_Data,type,solution_num,"axes",ax,"colour",colour_num,"backbone",0);
        end
        plot_backbone(Dyn_Data,type,solution_num,"axes",ax,"colour",0);
    case "amplitude"
        

        num_L_modes = length(L_modes);
        num_colours = size(get_plot_colours,1);
        for iMode = 1:num_L_modes
            Dyn_Data = Dyn_Data.validate_solution(solution_num,L_modes(iMode));
            colour_num = mod(iMode-1,num_colours)+1;
            plot_h_predicition(Dyn_Data,type,solution_num,"colour",colour_num,"backbone",1);
        end
end
end