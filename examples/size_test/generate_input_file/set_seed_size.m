clear

seed_size = 0.00307;
% 0.01    --> 7,000
% 0.005   --> 24,500
% 0.00307   --> 100,000
% 0.002   --> 300,000
% 0.00129 --> 1,000,000
% 0.001   --> 2,000,000

%-----------
SYSTEM_NAME = "mems_arch";
DIR_DELIMINATOR = "\";


system("python mesh_arch.py " + seed_size);
%pyrunfile("mesh_arch.py"); problems with abaqus package?


current_directory = pwd;
current_path = split(current_directory,DIR_DELIMINATOR);
previous_directory = join(current_path(1:(end-1)),DIR_DELIMINATOR);
geometry_path = previous_directory{1} + DIR_DELIMINATOR + "geometry" + DIR_DELIMINATOR + SYSTEM_NAME;

if isfolder(geometry_path)
   rmdir(geometry_path,"s")
end
mkdir(geometry_path)
copyfile(SYSTEM_NAME + ".inp",geometry_path)

%-----------
G_ID =fopen(SYSTEM_NAME + ".inp");
geometry = textscan(G_ID,'%s','delimiter','\n');
fclose(G_ID);
geometry = geometry{1,1};

node_position = read_abaqus_node_position(geometry);
mesh_data = get_mesh_data(geometry);
mesh_data = mesh_data{1};

num_nodes = size(node_position,1);
num_dof = num_nodes*mesh_data.dimension;

%-----------
% plot_fe_mesh(geometry)
str = regexprep(string(num_dof),'(?<!\.\d*)\d{1,3}(?=(\d{3})+\>)','$&,');
sprintf(str + " dofs")

%-----------
