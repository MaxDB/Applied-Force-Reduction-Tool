clear
close all
set_visualisation_level(3)

system_name = "shallow_arch_0";
Dyn_Data = initalise_dynamic_data(system_name);
Model = Dyn_Data.Dynamic_Model.Model;
%------------------
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
dof.position = [320,8.32,16]; %µm
dof.direction = 2;
Additional_Output.dof = dof;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);


ref_solution.name = "shallow_arch_14";
Dyn_Data_Ref = initalise_dynamic_data(ref_solution.name);
%--------------

ref_solution.sol_num = 2;
Ref_Sol = Dyn_Data_Ref.load_solution(ref_solution.sol_num);
ref_solution.orbit_subset = get_orbit_subset(Ref_Sol,[1.0279, 1.0341],3);

Force_Data = Ref_Sol.Force_Data;
Damping_Data = Ref_Sol.Damping_Data;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------

ref_solution.sol_num = 3;
Ref_Sol = Dyn_Data_Ref.load_solution(ref_solution.sol_num);
ref_solution.orbit_subset = get_orbit_subset(Ref_Sol,[1.0279, 1.0341],3);

Force_Data = Ref_Sol.Force_Data;
Damping_Data = Ref_Sol.Damping_Data;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------

ref_solution.sol_num = 4;
Ref_Sol = Dyn_Data_Ref.load_solution(ref_solution.sol_num);
ref_solution.orbit_subset = get_orbit_subset(Ref_Sol,[0.5062, 0.5269],3);

Force_Data = Ref_Sol.Force_Data;
Damping_Data = Ref_Sol.Damping_Data;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------

function orbit_subset = get_orbit_subset(Sol,span,interval)
freq = Sol.frequency;
include_orbit = (freq > span(1)) & (freq < span(2));
start_index = find(include_orbit,1);
end_index = find(include_orbit,1,"last");

orbit_subset = start_index:interval:end_index;

% test
% amp = Sol.energy;
% plot(freq(orbit_subset),amp(orbit_subset),"x");
end