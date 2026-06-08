clear
close all
set_visualisation_level(3)

system_name = "clamped_beam_0";
Dyn_Data = initalise_dynamic_data(system_name);
%------------------
Additional_Output.output = "physical displacement";
Additional_Output.type = "max";
Additional_Output.dof = 182;
Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);


Damping_Data.damping_type = "rayleigh";
Damping_Data.mass_factor = 0;
Damping_Data.stiffness_factor = 1.1e-3;

Force_Data.type = "point";
Force_Data.dof = 236;
Force_Data.continuation_variable = "frequency";
Force_Data.amplitude = 1;

%--------------
Dyn_Data_Ref = initalise_dynamic_data("clamped_beam_11001");
Sol = Dyn_Data_Ref.load_solution(3);
Force_Data.frequency = Sol.frequency;
%--------------

Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data);