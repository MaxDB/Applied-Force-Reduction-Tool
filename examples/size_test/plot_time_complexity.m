clear
data_name = ["none","perturbation"];
num_data = length(data_name);


load("data\size_data\size_data_none","Size_Data");
time = Size_Data.static_time;
time_all([1,3],:) = time;
num_dof = Size_Data.num_dof;

load("data\size_data\size_data_perturbation","Size_Data");
time = Size_Data.static_time;
time_all([2,4],:) = time;

figure
ax = axes();

box on
bar(log10(num_dof),time_all'/60,'stacked')
xlabel("log_{10}(degrees of freedom)")
ylabel("time (min)")

legend("$\mathcal W_{\mathcal R_1}$","$\mathcal V_{\mathcal R_1 : \mathcal R_{19}}$","$\mathcal W_{\mathcal R_2}$","$\mathcal V_{\mathcal R_2 : \mathcal R_{19}}$","Interpreter","latex")
