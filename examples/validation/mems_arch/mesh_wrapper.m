function num_dof = mesh_wrapper(seed_size)
%wrapper to get round live scripts not being able to call private functions
    num_dof = create_mesh(seed_size);
end