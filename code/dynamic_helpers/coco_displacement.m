function [data,disp_output] = coco_displacement(prob,data,u,disp_func)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
disp_output = disp_func(r); 
end