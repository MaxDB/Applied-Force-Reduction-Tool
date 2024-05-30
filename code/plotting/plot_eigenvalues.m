function plot_eigenvalues(Dyn_Data,solution_num,orbit_num)
NUM_THETA = 200;

if isstring(Dyn_Data)
    Dyn_Data = initalise_dynamic_data(Dyn_Data);
end


Rom = Dyn_Data.Dynamic_Model;
solution_name = Rom.data_path + "dynamic_sol_" + solution_num;
bd = coco_bd_read(solution_name);
solution_data  = coco_bd_col(bd, {"eigs"});
sol_labels = Dyn_Data.solution_labels{1,solution_num};
eigenvalues = solution_data(:,sol_labels(orbit_num));
ev_x = real(eigenvalues);
ev_y = imag(eigenvalues);

theta = linspace(0,2,NUM_THETA);
x = cospi(theta);
y = sinpi(theta);
figure
box on
hold on
plot(x,y,'k--')
plot(ev_x,ev_y,'x','MarkerSize',8,'LineWidth',1.5)
hold off
ylabel("Imaginary")
xlabel("Real")
end
