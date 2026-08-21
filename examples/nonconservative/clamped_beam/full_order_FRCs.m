clear
close all
set_visualisation_level(3)

system_name = "clamped_beam_0";
Dyn_Data = initalise_dynamic_data(system_name);
%------------------
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 248;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

%-----------------------------------------%
Force_Data.type = "point";
Force_Data.dof = 248;
Force_Data.continuation_variable = "frequency";

Damping_Data.damping_type = "rayleigh";
%-----------------------------------------%


ref_solution.name = "clamped_beam_11001";
%--------------
Force_Data.amplitude = 5;
target_damping = 0.3;
%-
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,0]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = 0;


ref_solution.sol_num = 1;
Dyn_Data_Ref = initalise_dynamic_data(ref_solution.name);
Ref_Sol = Dyn_Data_Ref.load_solution(ref_solution.sol_num);
ref_solution.orbit_subset = get_orbit_subset(Ref_Sol,[300,600],3);

Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);


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