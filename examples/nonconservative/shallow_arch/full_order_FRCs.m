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



Damping_Data.damping_type = "rayleigh";
damping_coeffs = get_shallow_arch_damping(Model);
Damping_Data.mass_factor = damping_coeffs(1);
Damping_Data.stiffness_factor = damping_coeffs(2);


Force_Data.type = "shape";
Force_Data.amplitude = 1;
kappa_1 = 0.03; %(µm/µs^2)
force_shape = @(kappa_2) get_shallow_arch_force(Model,kappa_1,kappa_2);

ref_solution.name = "shallow_arch_14";
%--------------

kappa_2 = 10;
Force_Data.shape = force_shape(kappa_2);
ref_solution.sol_num = 2;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------

kappa_2 = -10;
Force_Data.shape = force_shape(kappa_2);
ref_solution.sol_num = 3;
Dyn_Data = Dyn_Data.add_full_order_forced_response(Force_Data,Damping_Data,"solution",ref_solution);
%--------------


%--------------

function force_shape = get_shallow_arch_force(Model,kappa_1,kappa_2)
mass = Model.mass;
stiffness = Model.stiffness;

[evec,~] = eigs(stiffness,mass,4,"smallestabs");
bending_13 = evec(:,[1,4]);
force_shape = mass*bending_13*[kappa_1;kappa_2];
end

function damping_coeffs = get_shallow_arch_damping(Model)
mass = Model.mass;
stiffness = Model.stiffness;

[~,eval] = eigs(stiffness,mass,1,"smallestabs");

freq_1 = sqrt(eval);

damping_coeffs(1) = freq_1/500;
damping_coeffs(2) = 0;
end