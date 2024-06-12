function compare_validation(Dyn_Data,solution_num,L_modes)
figure
ax = gca;

num_L_modes = length(L_modes);
for iMode = 1:num_L_modes
    Dyn_Data = Dyn_Data.validate_solution(solution_num,L_modes(iMode));
    ax = plot_h_predicition(Dyn_Data,"energy",solution_num,ax);
end
plot_backbone(Dyn_Data,"energy",solution_num,"axes",ax);
end