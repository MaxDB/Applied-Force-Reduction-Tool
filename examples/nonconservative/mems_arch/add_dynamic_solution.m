clear

close all
set_visualisation_level(1)
set_logging_level(3)

system_name = "mems_arch_16";
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
Continuation_Opts.forward_steps = 0;
Continuation_Opts.backward_steps = 200;
Continuation_Opts.collocation_degree = 8;
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
[orbit,validation_orbit] = Dyn_Data_One_Mode.get_orbit(1,unstable_index(2),1);
[q,q_dot] = Dyn_Data_One_Mode.get_modal_validation_orbit(1,unstable_index(2));
[min_ke,min_index] = min(sum(q_dot.^2,1)); 
test_ic = q(:,min_index);

potential_ic = initial_condition_sweep(Dyn_Data.Dynamic_Model,2*pi/orbit.T,test_ic);

Continuation_Opts.initial_inc = 5e-1;
Continuation_Opts.max_inc = 5e-1;
Dyn_Data = Dyn_Data.add_backbone(1,"ic",potential_ic,"opts",Continuation_Opts);


%------------
orbit_id = 26;

figure;
ax = gca;
compare_orbits(["t","v-d-3"],"mems_arch_16",{2,orbit_id},"axes",ax,"colour",3);
compare_orbits(["t","q-d-5"],"mems_arch_16",{2,orbit_id},"axes",ax,"colour",0);

%------------
%for thesis
Model = Dyn_Data.Dynamic_Model.Model;
freq_1 = sqrt(Model.reduced_eigenvalues(1));
%---

Continuation_Opts.collocation_degree = 6;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.initial_discretisation_num = 40;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 500;
Continuation_Opts.frequency_range = [2.5,2.8]*1e6;

%-----------------------------------------%
Continuation_Opts.energy_limit_multiplier = 1;

Force_Data.type = "modal";
Force_Data.mode_number = 1;
% Force_Data.frequency = 2.72e6;
 Force_Data.frequency = 2.55e6;
  % Force_Data.frequency = 2.75e6;
Force_Data.amplitude = 20e3;
Force_Data.continuation_variable = "frequency";

%---
Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = freq_1/500;
Damping_Data.stiffness_factor = 0;
%---

Dyn_Data = Dyn_Data.add_forced_response(Force_Data,Damping_Data,"opts",Continuation_Opts);


