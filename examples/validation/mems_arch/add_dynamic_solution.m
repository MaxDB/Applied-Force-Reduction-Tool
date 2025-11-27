clear

close all
set_visualisation_level(1)
set_logging_level(3)

system_name = "mems_arch_1";
Dyn_Data = initalise_dynamic_data(system_name);
%-------------------------------------------------------------------------%
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [0,36e-3,10e-3];
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);
% --------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 5e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
switch system_name
    case "mems_arch_1"
        Continuation_Opts.initial_discretisation_num = 20;
    case "mems_arch_16"
        Continuation_Opts.initial_discretisation_num = 40; 
end
% -----------------------------------------%

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);

if system_name == "mems_arch_1", return, end

Dyn_Data_One_Mode = initalise_dynamic_data("mems_arch_1");
Dyn_Data_One_Mode = Dyn_Data_One_Mode.validate_solution(1,6);
Sol = Dyn_Data_One_Mode.load_solution(1,"validation");
unstable_index = find(Sol.h_stability>1.01,3);
[orbit,validation_orbit] = Dyn_Data_One_Mode.get_orbit(1,unstable_index(1),1);
[q,q_dot] = Dyn_Data_One_Mode.get_modal_validation_orbit(1,unstable_index(1));
[min_ke,min_index] = min(sum(q_dot.^2,1)); 
test_ic = q(:,min_index);

potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2*pi/orbit.T,test_ic);


Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);


%------------

