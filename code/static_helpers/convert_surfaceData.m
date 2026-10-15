function Static_Data = convert_surfaceData(Static_Data,Surface_Data)
Model = Static_Data.Model;

modal_force = Surface_Data.P;
modal_disp = Surface_Data.Q;
energy = Surface_Data.E;
energy = energy - Surface_Data.initialEnergy;

phy_disp = Surface_Data.eigenvectors*modal_disp;

r_disp = Model.reduced_eigenvectors'*Model.mass*phy_disp;
sep_ids = ones(size(energy));

Static_Data = Static_Data.update_data(r_disp,phy_disp,modal_force,energy,sep_ids,[]);

% ax = plot_static_data("displacement",Static_Data,"outputs",728);
% ax = ax{1};
% Rom = Reduced_System(Static_Data,"degree",[5,5]);
% hold(ax,"on")
% Rom.Physical_Displacement_Polynomial.plot_polynomial("outputs",728,"axes",ax)
% hold(ax,"off")
end