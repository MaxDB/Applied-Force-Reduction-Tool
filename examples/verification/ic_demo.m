clear
close all

spring_stiffness = [0,100,1000,1e4,1e10];
energy_limit = [0.18, 0.18, 0.18, 0.1, 0.01];

%--------- Software Settings ---------%
set_logging_level(1)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "ic_demo";
initial_modes = [1,2];
% initial_modes = 1;
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_loadcases = 10;
Static_Opts.maximum_loadcases = 20;
%------------------------------------------%

%--------- Static Verification Settings ---------%
Verification_Opts.verification_algorithm = "sep_to_edge";
Verification_Opts.maximum_iterations = 3;
Verification_Opts.maximum_interpolation_error = [1e-3,1e-2];
Verification_Opts.num_added_points = 1;
%----------------------------------------------%

Additional_Output.output = "physical displacement";
Additional_Output.type = "amplitude";
Additional_Output.dof = 362;

%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 1e-1;
Continuation_Opts.max_inc = 2e-1;
Continuation_Opts.min_inc = 1e-2;
Continuation_Opts.forward_steps = 100;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%

dynamic_system_name = system_name + "_" + join(string(initial_modes),"");
num_stiffness = length(spring_stiffness);

figure
ax = axes;

for iK = 1:num_stiffness
    stiffness_tag = string(iK);
    set_bc_stiffness(system_name,spring_stiffness(iK))
    Model = Dynamic_System(system_name,energy_limit(iK),initial_modes,Calibration_Opts,Static_Opts);

    Static_Data = Static_Dataset(Model,Verification_Opts);
    Static_Data.save_data;

    Dyn_Data = initalise_dynamic_data(dynamic_system_name);
    Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

    Continuation_Opts.inertial_compensation = 1;
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    ax = plot_backbone(Dyn_Data,"physical amplitude",1,"axes",ax,"tag",stiffness_tag + "_ic","colour",1);

    Continuation_Opts.inertial_compensation = 0;
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    ax = plot_backbone(Dyn_Data,"physical amplitude",2,"axes",ax,"tag",stiffness_tag,"colour",2);
end