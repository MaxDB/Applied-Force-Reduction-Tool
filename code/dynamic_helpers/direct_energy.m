function [data,E] = direct_energy(prob,data,u,Eom_Input)
pr = data.pr;
maps = pr.coll_seg.maps;

x = u(pr.xbp_idx);
xcn = reshape(maps.W*x, maps.x_shp);

num_modes = size(xcn,1)/2;

r = xcn(1:num_modes,:); %displacement
r_dot = xcn((num_modes+1):end,:);   %velocity

num_x = size(r,2);
kinetic_energy = zeros(1,num_x);
potential_energy = zeros(1,num_x);
for iX = 1:num_x
    kinetic_energy(1,iX) = 0.5*r_dot(:,iX)'*r_dot(:,iX);
    potential_energy(1,iX) =  Eom_Input.modal_potential(r(:,iX));
end

E = kinetic_energy+potential_energy;

if std(E)/mean(E) > 1e-2
    warning("energy not conserved")
end

E = max(E);
% data.energy = y;
% 
% y = max(T);
end