function C = get_rayleigh_damping_matrix(Damping_Data,Model)
    C = sparse(Damping_Data.mass_factor*Model.mass + Damping_Data.stiffness_factor*Model.stiffness);
end