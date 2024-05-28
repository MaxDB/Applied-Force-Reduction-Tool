function compare_solutions(system_name_1,system_name_2,sol_number_1,sol_number_2,type,validation)

if ~exist("validation","var")
    validation = 0;
end

Dyn_Data_1 = initalise_dynamic_data(system_name_1);
Dyn_Data_2 = initalise_dynamic_data(system_name_2);

switch validation
    case 0
        ax = plot_backbone(Dyn_Data_1,type,sol_number_1);
        plot_backbone(Dyn_Data_2,type,sol_number_2,ax);
    case 1
        ax = plot_h_predicition(Dyn_Data_1,type,sol_number_1);
        plot_backbone(Dyn_Data_2,type,sol_number_2,ax);
    case 2
        ax = plot_h_predicition(Dyn_Data_1,type,sol_number_1);
        plot_h_predicition(Dyn_Data_2,type,sol_number_2,ax);
end

end