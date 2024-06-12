function [data,phase] = coco_phase(prob,data,u,Applied_Force_Data)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

t = pr.coll_seg.mesh.tcn;

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
r_dot = xcn((num_modes+1):end,:);   %velocity

phase = 1;
end