function stress = solution_max_disp_stress(orbits,Rom,Additional_Output)
STRESS_TYPE = "von mises";
CALCULATION_TYPE = "extreme displacement";
STRESS_SECTION_POINT = 3;
num_dimensions = get_num_node_dimensions(Rom.Model);

Model = Rom.Model;
num_dofs = Model.num_dof;
num_modes = length(Model.reduced_modes);
node_map = Model.node_mapping;
num_nodes = max(node_map(:,1))/num_dimensions;


num_orbits = size(orbits,1);
if num_orbits == 1
    orbits = {orbits};
end


dof_bc = Additional_Output.dof;


dof = node_map(node_map(:,1) == dof_bc,2);
dof_node = 1 + (dof_bc - mod(dof_bc,6))/num_dimensions;

switch CALCULATION_TYPE
    case "extreme displacement"

        num_simulations = 2*num_orbits;
        physical_displacement = zeros(num_dofs,num_simulations);
        physical_velocity = zeros(num_dofs,num_simulations);
        initial_force = zeros(num_modes,num_simulations);

        for iOrbit = 1:num_orbits
            sim_index = 2*iOrbit - [1,0];

            orbit = orbits{iOrbit,1};
            r = orbit.xbp';


            r_disp = r(1:num_modes,:);
            x_disp = Rom.expand(r_disp);
            x_dof = x_disp(dof,:);

            [~,max_index] = max(x_dof);
            [~,min_index] = min(x_dof);

            physical_displacement(:,sim_index(1)) = x_disp(:,max_index);
            physical_displacement(:,sim_index(2)) = x_disp(:,min_index);

            r_disp_max = r(1:num_modes,max_index);
            r_disp_min = r(1:num_modes,min_index);

            r_vel_max = r((1:num_modes) + num_modes,max_index);
            r_vel_min = r((1:num_modes) + num_modes,min_index);

            initial_force(:,sim_index(1)) = Rom.Force_Polynomial.evaluate_polynomial(r_disp_max);
            initial_force(:,sim_index(2)) = Rom.Force_Polynomial.evaluate_polynomial(r_disp_min);

            physical_velocity(:,sim_index(1)) = Rom.expand_velocity(r_disp_max,r_vel_max);
            physical_velocity(:,sim_index(2)) = Rom.expand_velocity(r_disp_min,r_vel_min);

        end
    case "max displacement"

        num_simulations = num_orbits;
        physical_displacement = zeros(num_dofs,num_simulations);
        physical_velocity = zeros(num_dofs,num_simulations);
        initial_force = zeros(num_modes,num_simulations);

        for iOrbit = 1:num_orbits

            orbit = orbits{iOrbit,1};
            r = orbit.xbp';


            r_disp = r(1:num_modes,:);
            x_disp = Rom.expand(r_disp);
            x_dof = x_disp(dof,:);

            [~,max_index] = max(abs(x_dof));

            physical_displacement(:,iOrbit) = x_disp(:,max_index);

            r_disp_max = r(1:num_modes,max_index);

            r_vel_max = r((1:num_modes) + num_modes,max_index);

            initial_force(:,iOrbit) = Rom.Force_Polynomial.evaluate_polynomial(r_disp_max);

            physical_velocity(:,iOrbit) = Rom.expand_velocity(r_disp_max,r_vel_max);

        end
end

[stress_all,stress_labels,section_points] = Model.stress_simulation(physical_displacement,physical_velocity,initial_force);

stress = zeros(num_nodes,num_orbits);

switch CALCULATION_TYPE
    case "extreme displacement"
        for iOrbit = 1:num_orbits
            sim_index = 2*iOrbit - [1,0];

            section_point = section_points{1,sim_index(1)};
            stress_label = stress_labels{1,sim_index(1)};

            max_stress = stress_all{1,sim_index(1)}(:,:,section_point == STRESS_SECTION_POINT);
            min_stress = stress_all{1,sim_index(2)}(:,:,section_point == STRESS_SECTION_POINT);

            max_stress_tensor = construct_stress_tensor(max_stress,stress_label);
            min_stress_tensor = construct_stress_tensor(min_stress,stress_label);


            switch STRESS_TYPE
                case "von mises"
                    max_von_mises = get_von_mises_stress(max_stress_tensor);
                    min_von_mises = get_von_mises_stress(min_stress_tensor);
                    % [~,largest_stress_index] = max([max_von_mises(dof_node),min_von_mises(dof_node)]);
                    [~,largest_stress_index] = max([max(max_von_mises(dof_node)),max(min_von_mises(dof_node))]);
                    if largest_stress_index == 1
                        stress(:,iOrbit) = max_von_mises;
                    else
                        stress(:,iOrbit) = min_von_mises;
                    end
            end
        end
    case "max displacement"
        for iOrbit = 1:num_orbits


            section_point = section_points{1,iOrbit};
            stress_label = stress_labels{1,iOrbit};

            max_stress = stress_all{1,iOrbit}(:,:,section_point == STRESS_SECTION_POINT);
            max_stress_tensor = construct_stress_tensor(max_stress,stress_label);


            switch STRESS_TYPE
                case "von mises"
                    max_von_mises = get_von_mises_stress(max_stress_tensor);
                    stress(:,iOrbit) = max_von_mises;

            end
        end
end


end

function stress_tensor = construct_stress_tensor(stress,stress_label)
NUM_LINEAR_DIMENSIONS = 3;

num_nodes = size(stress,1);
stress_tensor = zeros(num_nodes,NUM_LINEAR_DIMENSIONS,NUM_LINEAR_DIMENSIONS);
num_components = size(stress,2);

for iComponent = 1:num_components
    label = strip(stress_label(iComponent),"left","S");
    label_array = convertStringsToChars(label);
    row_index = double(string(label_array(1)));
    col_index = double(string(label_array(2)));

    stress_tensor(:,row_index,col_index) = stress(:,iComponent);
    if row_index ~= col_index
        stress_tensor(:,col_index,row_index) = stress(:,iComponent);
    end

end

end

function von_mises_stress = get_von_mises_stress(stress)
direct_stress = 1/2*((stress(:,1,1) - stress(:,2,2)).^2 + (stress(:,2,2)-stress(:,3,3)).^2+(stress(:,3,3)-stress(:,1,1)).^2);
shear_stress = 3*(stress(:,2,3).^2+stress(:,3,1).^2 + stress(:,1,2).^2);
von_mises_stress = sqrt(direct_stress+shear_stress);
end