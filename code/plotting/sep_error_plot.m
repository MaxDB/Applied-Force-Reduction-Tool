function sep_error_plot(Plot_Data,Static_Data,Plot_Error,Disp_Error_Inputs)
PLOT_LEVEL = 4;
load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end
%------
fig = figure;
tiledlayout(2,2);

potential_ax = nexttile;
acceleration_ax = nexttile;
sep_ax = nexttile;
error_ax = nexttile;
%-----

disp = Plot_Data(1).displacement;
lambda = Plot_Data(1).lambda;
force_ratio = Plot_Data(1).force_ratio;
force = force_ratio.*lambda;

Rom_One = Plot_Data(1).Rom;
Rom_Two = Plot_Data(2).Rom;


%----
potential_one = Rom_One.Potential_Polynomial.evaluate_polynomial(disp);
potential_two = Rom_Two.Potential_Polynomial.evaluate_polynomial(disp);


set(fig, 'currentaxes', potential_ax);  
hold on
plot(lambda,potential_one)
plot(lambda,potential_two,"--")
hold off


title("n_f = " + Rom_One.Force_Polynomial.polynomial_degree + " : " + Rom_Two.Force_Polynomial.polynomial_degree)
box on
xlabel("\lambda")
ylabel("Potential")

%----
[r_ddot_one,r_ddot_two] = get_acceleration(disp,Rom_One,Rom_Two,Disp_Error_Inputs);


set(fig, 'currentaxes', acceleration_ax);  
hold on
plot(lambda,r_ddot_one)
plot(lambda,r_ddot_two,"--")
hold off


title("n_d = " + Rom_One.Physical_Displacement_Polynomial.polynomial_degree + " : " + Rom_Two.Physical_Displacement_Polynomial.polynomial_degree)
box on
xlabel("\lambda")
ylabel("Acceleration")


%----
num_modes = size(Static_Data,1);
force_points = Static_Data.restoring_force;

set(fig, 'currentaxes', sep_ax);  
hold on
switch num_modes
    case 1

    case 2
        plot(force_points(1,:),force_points(2,:),"x")
        plot(force(1,:),force(2,:))
    case 3
        plot3(force_points(1,:),force_points(2,:),force_points(3,:),"x")
        plot3(force(1,:),force(2,:),force(3,:))
        sep_ax.CameraPosition = [-1429.01208184779	-15572.6982084683	8848.89938869558];
end
hold off

box on
xlabel("f_1")
ylabel("f_2")
zlabel("f_3")


%----
tested_index = Plot_Data(1).tested_index;
tested_lambda = lambda(tested_index);
force_error = Plot_Error.force_error;
disp_error = Plot_Error.disp_error;

set(fig, 'currentaxes', error_ax);  
hold on
plot(tested_lambda ,force_error,".-")
plot(tested_lambda ,disp_error,".-")
plot(tested_lambda ,max(force_error,disp_error),"--")
plot(error_ax.XLim,[1,1],"k--")
hold off

legend("force","displacement")
box on
xlabel("\lambda")
ylabel("\epsilon")

%---
drawnow
end

function [r_ddot_one,r_ddot_two] = get_acceleration(disp,Rom_One,Rom_Two,Disp_Error_Inputs)
    beta_bar_one = Disp_Error_Inputs.beta_bar_one;
    beta_bar_two = Disp_Error_Inputs.beta_bar_two;
    input_order = Disp_Error_Inputs.input_order;
    Disp_Diff_One = Disp_Error_Inputs.Disp_Diff_Data_One;
    Disp_Diff_Two = Disp_Error_Inputs.Disp_Diff_Data_Two;

    scale_factor = Rom_One.Force_Polynomial.scaling_factor;
    shift_factor = Rom_One.Force_Polynomial.shifting_factor;
    r_transformed = scale_factor.*(disp + shift_factor);

    num_disp_coeffs = [size(beta_bar_one,1),size(beta_bar_two,1)];
    num_force_coeffs = [Rom_One.Force_Polynomial.num_element_coefficients,Rom_Two.Force_Polynomial.num_element_coefficients];
    num_coeffs = max(max(num_disp_coeffs),max(num_force_coeffs));

    num_x = size(r_transformed,2);
    num_modes = size(r_transformed,1);

    r_ddot_one = zeros(num_modes,num_x);
    r_ddot_two = zeros(num_modes,num_x);
    for iX = 1:num_x
        r_i = r_transformed(:,iX);

        r_power_products = ones(num_coeffs,1);
        for iMode = 1:num_modes
            r_power_products = r_power_products.*r_i(iMode).^input_order(:,iMode);
        end
        
        r_products_disp_one = r_power_products(1:num_disp_coeffs(1),:);
        r_products_disp_two = r_power_products(1:num_disp_coeffs(2),:);

        r_dr_products_disp_one = r_products_disp_one(Disp_Diff_One.diff_mapping{1,1}).*Disp_Diff_One.diff_scale_factor{1,1};
        r_dr_products_disp_two = r_products_disp_two(Disp_Diff_Two.diff_mapping{1,1}).*Disp_Diff_Two.diff_scale_factor{1,1};

        r_products_force_one = r_power_products(1:num_force_coeffs(1),:);
        r_products_force_two = r_power_products(1:num_force_coeffs(2),:);
    
        %----------
        force_one = Rom_One.Force_Polynomial.coefficients*r_products_force_one;
        force_two = Rom_Two.Force_Polynomial.coefficients*r_products_force_two;

        %---------
        inertia_one = r_dr_products_disp_one'*beta_bar_one*r_dr_products_disp_one;
        inertia_two = r_dr_products_disp_two'*beta_bar_two*r_dr_products_disp_two;

        %---------
        r_ddot_one(:,iX) = - inertia_one\force_one;
        r_ddot_two(:,iX) = - inertia_two\force_two;
    end
end