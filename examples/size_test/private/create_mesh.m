function num_dof = create_mesh(seed_size)
SYSTEM_NAME = "mems_arch";
DIR_DELIMINATOR = "\";

current_directory = pwd;
geometry_path = string(current_directory) + DIR_DELIMINATOR + "geometry" + DIR_DELIMINATOR + SYSTEM_NAME + DIR_DELIMINATOR;

reset_temp_directory()
cd temp
system("python " + '"' + geometry_path +  "mesh_arch.py" + '" ' + seed_size);
copyfile(SYSTEM_NAME + ".inp",geometry_path)
cd(current_directory)

%--
files_to_reset = {"matrices.mat","force_calibration.mat","mesh_data.mat"};
num_files = length(files_to_reset);
for iFile = 1:num_files
    file = geometry_path + files_to_reset{iFile};
    if isfile(file)
        delete(file)
    end
end

%--
G_ID =fopen(geometry_path + SYSTEM_NAME + ".inp");
geometry = textscan(G_ID,'%s','delimiter','\n');
fclose(G_ID);
geometry = geometry{1,1};

node_position = read_abaqus_node_position(geometry);
mesh_data = get_mesh_data(geometry);
mesh_data = mesh_data{1};

num_nodes = size(node_position,1);
num_dof = num_nodes*mesh_data.dimension;
end