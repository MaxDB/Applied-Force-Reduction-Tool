function [data,disp_max] = coco_displacement(prob,data,u,Theta_Poly,dof)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
phy_disp = Theta_Poly.evaluate_polynomial(r,dof); 
% disp_amplitude = (max(phy_disp) - min(phy_disp))/2;
disp_max = max(abs(phy_disp));
end