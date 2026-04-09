function final_force = rotate_force(starting_force,starting_displacement,final_displacement,Mesh_Data)

num_nodes = Mesh_Data.num_nodes;
num_dof = num_nodes*Mesh_Data.dimension;
initial_node_position = Mesh_Data.node_starting_position;

switch Mesh_Data.type
    case {"beam","shell"}
        if Mesh_Data.model_dimension == 3
            disp_dofs = 1:3;
            rotational_dofs = 4:6;
        end
end

nodal_dofs = ((1:Mesh_Data.dimension)+((0:(num_nodes-1))*Mesh_Data.dimension)')';


change_displacement = final_displacement - starting_displacement;

nodal_disp_all = change_displacement(nodal_dofs);
nodal_forces = starting_force(nodal_dofs);
rotated_nodal_forces = zeros(size(nodal_forces));
for iNode = 1:num_nodes
    nodal_rotation = nodal_disp_all(rotational_dofs,iNode);
    rotation_matrix = abaqus_rotation(nodal_rotation);

    force_disp = nodal_forces(disp_dofs,iNode);
    rotated_force_disp = rotation_matrix*force_disp;

    moments = nodal_forces(rotational_dofs,iNode);
    rotated_moments = sum(rotation_matrix*(diag(moments)),2);

    rotated_nodal_forces(:,iNode) = [rotated_force_disp;rotated_moments];

end

final_force = reshape(rotated_nodal_forces,num_dof,1);
end