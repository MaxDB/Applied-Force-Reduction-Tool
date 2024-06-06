function compare_solutions(type,varargin)


%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

dyn_data_names = cell(0,1);
solution_index = cell(0,1);
system_counter = 0;
for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "validation"
            validation = keyword_values{arg_counter};
        otherwise
            system_counter = system_counter + 1;
            dyn_data_names{1,system_counter} = keyword_args{arg_counter};
            solution_index{1,system_counter} = keyword_values{arg_counter};
    end
end
%-------------------------------------------------------------------------%
num_solutions = system_counter;
if ~exist("validation","var")
    validation = zeros(1,num_solutions);
elseif isscalar(validation)
    validation = validation * ones(1,num_solutions);
end

colour_numbers = 1:num_solutions;

ax = [];
for iSol = 1:num_solutions
    data_name = dyn_data_names{1,iSol};
    Dyn_Data = initalise_dynamic_data(data_name);

    system_sols = solution_index{1,iSol};
    num_system_sols = length(system_sols);
    for jSol = 1:num_system_sols
        if validation(num_solutions)
            ax = plot_h_predicition(Dyn_Data,type,system_sols(jSol),"axes",ax,"colour",colour_numbers(iSol));
        else
            ax = plot_backbone(Dyn_Data,type,system_sols(jSol),"axes",ax,"colour",colour_numbers(iSol));
        end
    end
end

end