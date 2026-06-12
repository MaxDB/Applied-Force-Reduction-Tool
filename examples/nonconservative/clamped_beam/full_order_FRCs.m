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


Damping_Data.damping_type = "rayleigh";
Force_Data.type = "point";
Force_Data.dof = 248;
Force_Data.continuation_variable = "frequency";

ref_solution.name = "clamped_beam_11001";
%--------------
Force_Data.amplitude = 0.05;
target_damping = 0.01;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);


ref_solution.sol_num = 2; %43
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------
Force_Data.amplitude = 1;
target_damping = 0.2;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);

ref_solution.sol_num = 3;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------
Force_Data.amplitude = 5;
target_damping = 1;
damping_coeffs = get_rayleigh_coeffs(Dyn_Data.Dynamic_Model.Model,target_damping,[1,2]);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);

ref_solution.sol_num = 4;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------


%--------------

