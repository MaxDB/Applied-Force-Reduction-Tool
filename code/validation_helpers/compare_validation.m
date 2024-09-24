function compare_validation(Dyn_Data,type,solution_num,h_modes)
ax = [];

if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end

if isstring(h_modes) && h_modes == "all"
    h_modes = 1:Dyn_Data.Dynamic_Model.Model.Static_Options.num_validation_modes; 
end
r_modes = Dyn_Data.Dynamic_Model.Model.reduced_modes;
L_modes = setdiff(h_modes,r_modes);

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