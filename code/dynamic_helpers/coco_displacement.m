function [data,disp_amplitude] = coco_displacement(prob,data,u,Theta_Poly,dof,r_evec)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
theta = Theta_Poly.evaluate_polynomial(r,dof);
phy_disp = r_evec(dof,:)*r + theta; 
disp_amplitude = (max(phy_disp) - min(phy_disp))/2;
end