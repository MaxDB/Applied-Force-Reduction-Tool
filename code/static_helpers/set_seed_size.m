function num_dof = set_seed_size(system_name,seed_size)
geometry_path = "geometry\" + system_name + "\";
cae_path = "\" + geometry_path + system_name + ".cae";
cae_path = replace(cae_path,"\","\\");

python_path = get_project_path + "\code\abaqus_scripts\set_seed_size.py";
python_path = replace(python_path,"\","/");

cmd_args = 'python "' + python_path + '" "' + cae_path + '" ' + seed_size;

reset_temp_directory()
current_dir = pwd;
cd("temp\")
[cmd_state,console_out] = system(cmd_args);
cd(current_dir)

reset_data(system_name)
geometry = geometry_path + system_name + ".inp";
movefile("temp\remesh.inp",geometry)

node_position = read_abaqus_node_position(geometry);
mesh_data = get_mesh_data(geometry);
mesh_data = mesh_data{1};

num_nodes = size(node_position,1);
num_dof = num_nodes*mesh_data.dimension;

end

