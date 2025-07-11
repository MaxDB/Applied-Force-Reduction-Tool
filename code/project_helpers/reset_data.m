function reset_data(system_name)
delete_static_data(system_name);
delete_cache(system_name,"force")
delete_cache(system_name,"matrices")
delete_cache(system_name,"mesh_data")
end