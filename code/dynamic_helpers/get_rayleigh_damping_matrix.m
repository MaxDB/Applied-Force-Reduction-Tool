function C = get_rayleigh_damping_matrix(Damping_Data,Model)
    C = Damping_Data.alpha*Model.mass + Damping_Data.beta*Model.stiffness;
end