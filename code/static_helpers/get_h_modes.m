function [h_evec,h_eval] = get_h_modes(Static_Data,h_modes)
Validation_Data = Static_Data.Dynamic_Validation_Data;
if ~isempty(Validation_Data)
    if Validation_Data.largest_h_mode >= max(h_modes)
        h_evec = Validation_Data.h_eigenvectors;
        h_eval = Validation_Data.h_eigenvalues;
        return
    end
end
Model = Static_Data.Model;
mass = Model.mass;
stiffness = Model.stiffness;

[h_evec,h_eval] = eigs(stiffness,mass,max(h_modes),"smallestabs");
h_eval = diag(h_eval);

r_modes = Model.reduced_modes;
h_eval(r_modes) = [];
h_evec(:,r_modes) = [];
end