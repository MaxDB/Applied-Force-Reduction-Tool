function num_dimensions = get_num_node_dimensions(Model)
mesh_data_path = "geometry\" + Model.system_name + "\mesh_data";
load(mesh_data_path,"mesh_data");

num_dimensions = mesh_data{1}.dimension;
end