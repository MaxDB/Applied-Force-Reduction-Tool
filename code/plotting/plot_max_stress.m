function plot_max_stress(Dyn_Data,solution_num,orbit_num,dof)
NUM_DIMENSIONS = 6;

Rom = Dyn_Data.Dynamic_Model;
Model = Rom.Model;
num_dof_bc = max(Model.node_mapping(:,1));
num_nodes = num_dof_bc/NUM_DIMENSIONS;

if nargin == 3
    % node_map = Model.node_mapping;
    % dof = node_map(node_map(:,1) == Dyn_Data.Additional_Output.dof,2);
    dof = Dyn_Data.Additional_Output.dof;
end

orbit = Dyn_Data.get_orbit(solution_num,orbit_num);

z = orbit.xbp';
num_modes = size(z,1)/2;
r = z(1:num_modes,:);
x = Rom.expand(r,"full");
x_dof = x(dof,:);
[~,max_index] = max(abs(x_dof));
x_max = x(:,max_index);

node_disp = reshape(x_max,NUM_DIMENSIONS,num_nodes)';
colour_map = Dyn_Data.max_displacement_stress{1,solution_num}(:,orbit_num);

plot_fe_mesh(Model,node_disp,colour_map)
end