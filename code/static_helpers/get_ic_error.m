function ic_error = get_ic_error(Static_Data)
    Rom = Reduced_System(Static_Data);
    mass = Rom.Model.mass;

    Disp_Poly = Rom.Physical_Displacement_Polynomial;
    Disp_dr_Poly = Disp_Poly.differentiate_polynomial;

    r = Static_Data.get_dataset_values("reduced_displacement");
    num_points = size(r,2);
    num_modes = size(r,1);
    I_r = eye(num_modes);
    ic_error = 0;
    for iPoint = 1:num_points
        r_i = r(:,iPoint);
        x_dr = Disp_dr_Poly.evaluate_polynomial(r_i);
        inertia_term = x_dr'*mass*x_dr;
        error = norm(inertia_term - I_r);
        ic_error = max(ic_error,error);
    end
end