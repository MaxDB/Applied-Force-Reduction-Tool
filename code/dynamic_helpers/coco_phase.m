function [data,phase] = coco_phase(prob,data,u,Applied_Force_Data)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

t = pr.coll_seg.mesh.tcn';
% t = pr.coll_seg.mesh.tbp(maps.tbp_idx);

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
% r_dot = xcn((num_modes+1):end,:);   %velocity

num_points = 10;
t_lin = linspace(t(1),1,num_points+1);
t_lin(end) = [];

L = length(t_lin);   %discrete length

mode_num = 1;
r_lin = interp1(t,r(mode_num,:),t_lin);
% r_signal = r_lin(1:(end-1)); %!!!!!!
r_signal = r_lin;
X = fft(r_signal);
X = X/(L-1);

alpha = 2*real(X(1,2));
beta = -2*imag(X(1,2));

phi = atan2(-beta,alpha);
if phi > 0
    phi = phi - 2*pi;
end
phase = phi + pi/2;
end
