function plot_eigenvalues(Dyn_Data,solution_num,orbit_num,varargin)
NUM_THETA = 200;
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
validated = 0;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "axes"
            ax = keyword_values{arg_counter};
        case {"validation"}
            validated = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%
if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end

if isempty(ax)
    figure
    ax = axes;
    box(ax,"on")
end


Rom = Dyn_Data.Dynamic_Model;
Sol = Dyn_Data.load_solution(solution_num);
if validated
    [~,validated_orbit] = Dyn_Data.get_orbit(solution_num,orbit_num,validated);
    validated_eigenvalues = validated_orbit.evals;
    validated_ev_x = real(validated_eigenvalues);
    validated_ev_y = imag(validated_eigenvalues);
end
solution_name = Rom.data_path + "dynamic_sol_" + solution_num;
bd = coco_bd_read(solution_name);
solution_data  = coco_bd_col(bd, {"eigs"});
orbit_label = Sol.orbit_labels(orbit_num);
eigenvalues = solution_data(:,orbit_label);
ev_x = real(eigenvalues);
ev_y = imag(eigenvalues);

theta = linspace(0,2,NUM_THETA);
x = cospi(theta);
y = sinpi(theta);


hold(ax,"on")
plot(ax,x,y,'k--')
plot(ax,ev_x,ev_y,'x','MarkerSize',8,'LineWidth',1.5)
if validated
    plot(ax,validated_ev_x,validated_ev_y,'x','MarkerSize',8,'LineWidth',1.5)
end
hold(ax,"on")
ylabel(ax,"Imaginary")
xlabel(ax,"Real")
end
