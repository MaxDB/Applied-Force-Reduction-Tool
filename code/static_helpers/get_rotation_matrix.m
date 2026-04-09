function rotation_matrix = get_rotation_matrix(angles)
num_dims = length(angles);

sin_angles = sin(angles);
cos_angles = cos(angles);

switch num_dims
    case 1
        rotation_matrix = [cos_angles,sin_angles;-sin_angles,cos_angles]; 
    case 3
        rotation_matrix = []

end

end