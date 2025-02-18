clear
close all

x_spring_stiffness = 1e3;
y_spring_stiffness = [1.78e3,2.86e3,6.08e3,1e8]; %0.05 -- 0.21
% y_spring_stiffness = 1e3;
% x_spring_stiffness = [1e4,3e4,5e4,7e4];
energy_limit = 2.5;

%--------- Software Settings ---------%
set_logging_level(1)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "ic_demo";
initial_modes = [1];
% initial_modes = 1;
%-----------------------------------%

%--------- Calibration Settings ---------%
Calibration_Opts.calibration_scale_factor = 2;
Calibration_Opts.Static_Opts.num_loadcases = 20;
%----------------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 4; %be careful!
Static_Opts.num_loadcases = 15;
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
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collation_degree = 6;
%-----------------------------------------%

dynamic_system_name = system_name + "_" + join(string(initial_modes),"");
num_stiffness = length(y_spring_stiffness);

fig_bb = figure;
ax_bb = axes(fig_bb);

fig_x = figure;
ax_x = axes(fig_x);

ic_error = zeros(1,num_stiffness);
reduced_modeshape = zeros(0,num_stiffness);
most_coupled_modeshape = zeros(0,num_stiffness);
most_coupled_modes = zeros(10,num_stiffness);
for iK = 1:num_stiffness
    stiffness_tag = string(iK);
    set_bc_stiffness(system_name,x_spring_stiffness,y_spring_stiffness(iK))
    Model = Dynamic_System(system_name,energy_limit,initial_modes,Calibration_Opts,Static_Opts);

    Static_Data = Static_Dataset(Model,Verification_Opts);
    ic_error(iK) = get_ic_error(Static_Data); 
    Static_Data.save_data;

    Dyn_Data = initalise_dynamic_data(dynamic_system_name);
    Dyn_Data = Dyn_Data.add_additional_output(Additional_Output);

    Continuation_Opts.inertial_compensation = 1;
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    ax_bb = plot_backbone(Dyn_Data,"physical amplitude",1,"axes",ax_bb,"tag",stiffness_tag + "_ic","colour",1);

    Continuation_Opts.inertial_compensation = 0;
    Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    ax_bb = plot_backbone(Dyn_Data,"physical amplitude",2,"axes",ax_bb,"tag",stiffness_tag,"colour",2);

    %plot quasi-static coupling
    Rom = Dyn_Data.Dynamic_Model;
    num_dof = Rom.Model.num_dof;
    [ranked_mode_list,ranked_modeshapes,Modal_Disp_Poly] = rank_condensed_modes(Rom);
    reduced_modeshape(1:num_dof,iK) = ranked_modeshapes(:,1);
    mode_sign = sign(reduced_modeshape(2,iK));
    reduced_modeshape(:,iK) = reduced_modeshape(:,iK)*mode_sign;

    most_coupled_modeshape(1:num_dof,iK) = ranked_modeshapes(:,2)*mode_sign;
    most_coupled_modes(:,iK) = ranked_mode_list(1:10,1);

    Modal_Disp_Poly = -1*Modal_Disp_Poly;
    Modal_Disp_Poly.plot_polynomial(ax_x,ranked_mode_list(2));
end