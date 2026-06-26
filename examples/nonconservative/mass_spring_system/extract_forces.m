clear
close all
comp_pair = 3;
%------

system_names = ["mass_spring_roller_0","mass_spring_roller_12"];
sol_nums = [2,1] + (comp_pair-1);
plot_lines = ["inertial","damping","elastic","external"];
% plot_lines = ["damping","external"];

num_systems = length(system_names);

force_fig = figure;
tiles = tiledlayout(force_fig,3,num_systems);

energy_fig = figure;
energy_ax = axes(energy_fig);
energy_plot_style = ["-","--"];
box(energy_ax,"on")


for iSystem = 1:num_systems
    system_name = system_names(iSystem);
    sol_num = sol_nums(iSystem);
    Dyn_Data = initalise_dynamic_data(system_name);
    Sol = Dyn_Data.load_solution(sol_num);
    orbit_id = Dyn_Data.get_special_point(sol_num,"RES");
    Orbit = Dyn_Data.get_orbit(sol_num,orbit_id);


    %----
    Forcing_Terms = separate_forcing_terms(Dyn_Data,sol_num,orbit_id);
    t = Orbit.tbp';
    %--

    f_inertia = Forcing_Terms.inertia_force + Forcing_Terms.conv_force;
    f_damping = Forcing_Terms.damping_force;
    f_restoring = Forcing_Terms.restoring_force;
    f_applied = -Forcing_Terms.applied_force;


    for iMode = 1:3
        tile_id = 2*(iMode-1) + iSystem;
        ax = nexttile(tiles,tile_id);
        hold(ax,"on")
        if any(f_inertia(iMode,:) ~= 0) && ismember("inertial",plot_lines)
            plot(ax,t,f_inertia(iMode,:),"Color",get_plot_colours(1))
        end
        if any(f_damping(iMode,:) ~= 0) && ismember("damping",plot_lines)
            plot(ax,t,f_damping(iMode,:),"Color",get_plot_colours(2))
        end
        if any(f_restoring(iMode,:) ~= 0) && ismember("elastic",plot_lines)
            plot(ax,t,f_restoring(iMode,:),"Color",get_plot_colours(3))
        end
        if any(f_applied(iMode,:) ~= 0) && ismember("external",plot_lines)
            plot(ax,t,f_applied(iMode,:),"Color",get_plot_colours(4))
        end
        hold(ax,"off")

        box(ax,"on")
        xlabel(ax,"Time (s)")
        ylabel(ax,"f_{q_" + iMode + "}")

        xlim(ax,[min(t),max(t)])

    end
    
    hold(energy_ax,"on")
    plot(energy_ax,t,Forcing_Terms.inertia_work_done,"LineStyle",energy_plot_style(iSystem),"Color",get_plot_colours(1))
    plot(energy_ax,t,Forcing_Terms.damp_work_done,"LineStyle",energy_plot_style(iSystem),"Color",get_plot_colours(2))
    plot(energy_ax,t,Forcing_Terms.elastic_work_done,"LineStyle",energy_plot_style(iSystem),"Color",get_plot_colours(3))
    plot(energy_ax,t,Forcing_Terms.ext_work_done,"LineStyle",energy_plot_style(iSystem),"Color",get_plot_colours(4))
    % plot(energy_ax,t,Forcing_Terms.ext_work_done + Forcing_Terms.damp_work_done,"LineStyle",energy_plot_style(iSystem))
    hold(energy_ax,"off")
    
    xlim(energy_ax,[min(t),max(t)])
end

hold(ax,"on")
if ismember("inertial",plot_lines)
    plot(ax,0,0,"Color",get_plot_colours(1))
end
if ismember("damping",plot_lines)
plot(ax,0,0,"Color",get_plot_colours(2))
end
if ismember("elastic",plot_lines)
plot(ax,0,0,"Color",get_plot_colours(3))
end
if ismember("external",plot_lines)
plot(ax,0,0,"Color",get_plot_colours(4))
end
hold(ax,"off")
legend(ax,plot_lines)


hold(energy_ax,"on")
lines(1) = plot(energy_ax,0,0,"Color",get_plot_colours(1));
lines(2) = plot(energy_ax,0,0,"Color",get_plot_colours(2));
lines(3) = plot(energy_ax,0,0,"Color",get_plot_colours(3));
lines(4) = plot(energy_ax,0,0,"Color",get_plot_colours(4));
lines(5) = plot(energy_ax,0,0,"Color",get_plot_colours("grey"));
lines(6) = plot(energy_ax,0,0,"--","Color",get_plot_colours("grey"));
hold(energy_ax,"off")


legend(lines,"Inertial","Damping","Elastic","External","FOM","ROM")

function Forcing_Terms = separate_forcing_terms(Dyn_Data,sol_num,orbit_id)
Sol = Dyn_Data.load_solution(sol_num);
Orbit = Dyn_Data.get_orbit(sol_num,orbit_id);


Rom = Dyn_Data.Dynamic_Model;
Model = Rom.Model;
[evec,eval] = eig(Model.stiffness,Model.mass);


Sol_Type = Dyn_Data.solution_types{sol_num};

Force_Data = [];
Damping_Data = [];
if Sol_Type.orbit_type == "forced"
    Force_Data = Sol_Type.Force_Data;
    Force_Data.frequency = 2*pi/Orbit.T;
    Damping_Data = Sol_Type.Damping_Data;
    Nonconservative_Input = Forced_Solution.get_nonconservative_input(Force_Data,Damping_Data,Rom);
end

eom = Rom.get_equation_of_motion("forcing",Force_Data,"damping",Damping_Data);


z = Orbit.xbp';
t = Orbit.tbp';
num_time_points = size(t,2);

z_dot = eom(t,z);

num_reduced_dofs = size(z,1)/2;
disp_index = 1:num_reduced_dofs;
vel_index = disp_index + num_reduced_dofs;

r = z(disp_index,:);
r_dot = z(vel_index,:);
r_ddot = z_dot(vel_index,:);

switch Sol_Type.model_type
    case "fom"
        model_dofs = Model.num_dof;
    case "rom"
        model_dofs = length(Model.reduced_modes);
end
inertia_force = zeros(model_dofs,num_time_points);
conv_force = zeros(model_dofs,num_time_points);
restoring_force = zeros(model_dofs,num_time_points);
damping_force = zeros(model_dofs,num_time_points);
applied_force = zeros(model_dofs,num_time_points);

inertia_work_done = zeros(1,num_time_points);
ext_work_done = zeros(1,num_time_points);
elastic_work_done = zeros(1,num_time_points);
damping_work_done = zeros(1,num_time_points);

switch Sol_Type.model_type
    case "fom"
        Analytic_Eom = load_analytic_system("geometry\" + Model.system_name+ "\" + Model.system_name);
        %--
        coord_transform = Analytic_Eom.linear_mass*evec;
        x = coord_transform*r;
        x_dot = coord_transform*r_dot;
        x_ddot = coord_transform*r_ddot;
        
        %--
        inertia_force = Analytic_Eom.linear_mass*x_ddot;
        %--
        restoring_force = Analytic_Eom.linear_stiffness*x;
        for iTime = 1:num_time_points
           restoring_force(:,iTime) = restoring_force(:,iTime) + Analytic_Eom.nonlinear_restoring_force(x(:,iTime));
        end
        %--
        damping_force = Nonconservative_Input.damping*x_dot;
        %--
        sin_term = sin(Force_Data.frequency*t);
        applied_force = Nonconservative_Input.amplitude*Nonconservative_Input.amplitude_shape.*sin_term;
        %--
        
        inertia_work = @(x)work_done(x,t,inertia_force,x_dot);
        damp_work = @(x)work_done(x,t,damping_force,x_dot);
        elastic_work = @(x)work_done(x,t,restoring_force,x_dot);
        ext_work = @(x)work_done(x,t,-applied_force,x_dot);
        t_start = t(1);
        total_inertia_work = 0;
        total_damping_work = 0;
        total_elastic_work = 0;
        total_ext_work = 0;
        for iTime = 1:num_time_points
            t_end = t(iTime);
            inertia_work_inc = integral(inertia_work,t_start,t_end);
            damping_work_inc = integral(damp_work,t_start,t_end);
            elastic_work_inc = integral(elastic_work,t_start,t_end);
            ext_work_inc = integral(ext_work,t_start,t_end);

            total_inertia_work = total_inertia_work + inertia_work_inc;
            total_damping_work = total_damping_work + damping_work_inc;
            total_elastic_work = total_elastic_work + elastic_work_inc;
            total_ext_work = total_ext_work + ext_work_inc;
            
            inertia_work_done(iTime) = total_inertia_work;
            damping_work_done(iTime) = total_damping_work;
            elastic_work_done(iTime) = total_elastic_work;
            ext_work_done(iTime) = total_ext_work;

            t_start = t_end;
        end
        


    case "rom"
        f_r = Rom.Force_Polynomial;
        x_dr = Rom.Physical_Displacement_Polynomial.differentiate_polynomial;
    
        sin_term = sin(Force_Data.frequency*t);
        phy_applied_force = Nonconservative_Input.amplitude*Nonconservative_Input.amplitude_shape.*sin_term;

        Eom_Data = Rom.get_solver_inputs("coco_backbone");
        Disp_Data = Eom_Data.Disp_Data;
        disp_coeffs = Rom.Physical_Displacement_Polynomial.coefficients;
        for iTime = 1:num_time_points
            r_i = r(:,iTime);
            r_dot_i = r_dot(:,iTime);
            r_ddot_i = r_ddot(:,iTime);

            x_dr_i = x_dr.evaluate_polynomial(r_i);
            f_r_i = f_r.evaluate_polynomial(r_i);

            inertia_force(:,iTime) = x_dr_i'*Model.mass*x_dr_i*r_ddot_i;
            %--
            r_power_products = ones(size(disp_coeffs,2),1);
            for iMode = 1:model_dofs
                r_power_products = r_power_products.*r_i(iMode).^Eom_Data.input_order(:,iMode);
            end
            
            r_dr2_products_disp = r_power_products(Disp_Data.diff_mapping{1,2}).*Disp_Data.diff_scale_factor{1,2};
            x_dr2_i = tensorprod(disp_coeffs,r_dr2_products_disp,2,1);
            
            conv_force(:,iTime) = x_dr_i'*Model.mass*tensorprod(x_dr2_i,r_dot_i,3,1)*r_dot_i;
            %--
            damping_force(:,iTime) = x_dr_i'*Nonconservative_Input.damping*x_dr_i*r_dot_i;
            %--
            restoring_force(:,iTime) = x_dr_i'*Model.mass*Model.reduced_eigenvectors*f_r_i;
            %--
            applied_force(:,iTime) = x_dr_i'*phy_applied_force(:,iTime);
        end
        %--
        x_dot = Rom.expand_velocity(r,r_dot);
        phy_transform = Model.mass*Model.reduced_eigenvectors;
        inertia_work = @(x)work_done(x,t,phy_transform*inertia_force,x_dot);
        damp_work = @(x)work_done(x,t,phy_transform*damping_force,x_dot);
        elastic_work = @(x)work_done(x,t,phy_transform*restoring_force,x_dot);
        ext_work = @(x)work_done(x,t,-phy_transform*applied_force,x_dot);

      
        t_start = t(1);
        total_inertia_work = 0;
        total_damping_work = 0;
        total_elastic_work = 0;
        total_ext_work = 0;
        for iTime = 1:num_time_points
            t_end = t(iTime);
            inertia_work_inc = integral(inertia_work,t_start,t_end);
            damping_work_inc = integral(damp_work,t_start,t_end);
            elastic_work_inc = integral(elastic_work,t_start,t_end);
            ext_work_inc = integral(ext_work,t_start,t_end);

            total_inertia_work = total_inertia_work + inertia_work_inc;
            total_damping_work = total_damping_work + damping_work_inc;
            total_elastic_work = total_elastic_work + elastic_work_inc;
            total_ext_work = total_ext_work + ext_work_inc;

            inertia_work_done(iTime) = total_inertia_work;
            damping_work_done(iTime) = total_damping_work;
            elastic_work_done(iTime) = total_elastic_work;
            ext_work_done(iTime) = total_ext_work;

            t_start = t_end;
        end


        %--
        
       



        phy_transform = Model.mass*Model.reduced_eigenvectors;
        inertia_force = phy_transform*inertia_force;
        conv_force = phy_transform*conv_force;
        damping_force = phy_transform*damping_force;
        restoring_force = phy_transform*restoring_force;
        applied_force = phy_transform*applied_force;
end




Forcing_Terms.inertia_force = evec'*inertia_force;
Forcing_Terms.conv_force = evec'*conv_force;
Forcing_Terms.restoring_force = evec'*restoring_force;
Forcing_Terms.damping_force = evec'*damping_force;
Forcing_Terms.applied_force = evec'*applied_force;


Forcing_Terms.inertia_work_done = inertia_work_done;
Forcing_Terms.damp_work_done = damping_work_done;
Forcing_Terms.elastic_work_done = elastic_work_done;
Forcing_Terms.ext_work_done = ext_work_done;

% Forcing_Terms.inertia_force = inertia_force;
% Forcing_Terms.conv_force = conv_force;
% Forcing_Terms.restoring_force = restoring_force;
% Forcing_Terms.damping_force = damping_force;
% Forcing_Terms.applied_force = applied_force;
end


function y = work_done(t,t_all,f_all,vel_all)
num_time_points = length(t);
num_dof = size(f_all,1);
f = zeros(num_dof,num_time_points);
vel = f;
for iDof = 1:num_dof
    f(iDof,:) = interp1(t_all,f_all(iDof,:),t,"spline");
    vel(iDof,:) = interp1(t_all,vel_all(iDof,:),t,"spline");
end
y = dot(f,vel);
end