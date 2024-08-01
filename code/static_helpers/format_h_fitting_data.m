function [h_input,h_output] = format_h_fitting_data(type,Static_Data,r,perturbation_displacement,base_output)

Model = Static_Data.Model;
r_evec = Model.reduced_eigenvectors;
L_evec = Static_Data.Dynamic_Validation_Data.current_L_eigenvectors;
h_evec = [r_evec,L_evec];

num_r_modes = size(r_evec,2);
num_h_modes = size(h_evec,2);
num_dofs = size(r_evec,1);

h_transform = (Model.mass*h_evec)';

num_loadcases = size(r,2);
loadcase_length = (num_h_modes+1);
h_input = zeros(num_r_modes + num_h_modes,num_loadcases*loadcase_length );


switch type
    case "displacement"
        h_output = zeros(num_dofs,num_loadcases*loadcase_length );
        for iLoad = 1:num_loadcases
            loadcase_span = ((0:(loadcase_length-1))*num_loadcases) + iLoad;

            perturbation_i = perturbation_displacement(:,:,iLoad);
            h_i = h_transform*perturbation_i;

            r_i = r(:,iLoad);

            h_input(1:num_r_modes,loadcase_span) = repmat(r_i,1,loadcase_length);
            h_input((num_r_modes+1):end,loadcase_span) = [zeros(num_h_modes,1),h_i];

            h_output(:,loadcase_span) = [zeros(num_dofs,1),perturbation_i] + base_output(:,iLoad);
        end
    case "force"
        h_output = zeros(num_r_modes + num_h_modes,num_loadcases*loadcase_length );
        h_force = Static_Data.perturbation_scale_factor*eye(num_h_modes);
        for iLoad = 1:num_loadcases
            loadcase_span = ((0:(loadcase_length-1))*num_loadcases) + iLoad;

            perturbation_i = perturbation_displacement(:,:,iLoad);
            h_i = h_transform*perturbation_i;

            r_i = r(:,iLoad);

            h_input(1:num_r_modes,loadcase_span) = repmat(r_i,1,loadcase_length);
            h_input((num_r_modes+1):end,loadcase_span) = [zeros(num_h_modes,1),h_i];
            
            h_output(1:num_r_modes,loadcase_span) = repmat(base_output(:,iLoad),1,loadcase_length);
            h_output((num_r_modes+1):end,loadcase_span) = [zeros(num_h_modes,1),h_force];
        end
end
end