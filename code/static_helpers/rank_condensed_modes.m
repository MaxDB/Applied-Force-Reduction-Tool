function [ranked_mode_list,ranked_modeshapes,Modal_Disp_Poly] = rank_condensed_modes(Rom)
    Model = Rom.Model;
    mass = Model.mass;
    stiffness = Model.stiffness;
    [evec_all,~] = eig(full(stiffness),full(mass),"vector");
    modal_transform = evec_all'*mass;
    Phy_Disp_Poly = Rom.Physical_Displacement_Polynomial;
    Modal_Disp_Poly = modal_transform*Phy_Disp_Poly;

    r_lim = Modal_Disp_Poly.input_limit;
    r = linspace(r_lim(1),r_lim(2),100);

    modal_disp = Modal_Disp_Poly.evaluate_polynomial(r);
    max_abs_disp = max(abs(modal_disp),[],2);
    [~,ranked_mode_list] = sort(max_abs_disp,"descend");
    ranked_modeshapes = evec_all(:,ranked_mode_list);
end