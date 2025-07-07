clear
num_evals = 6;
num_steps = 2;

max_seed_size = 0.01; %(7,000 dofs)
min_seed_size = 0.00129; %(1,000,000 dofs)

%-----------
seed_sizes = flip(logspace(log10(min_seed_size),log10(max_seed_size),num_steps));
SYSTEM_NAME = "mems_arch";
DIR_DELIMINATOR = "\";

num_dofs = zeros(1,num_steps);
natural_frequency = zeros(num_evals,num_steps);
for iStep = 1:num_steps
    seed_size = seed_sizes(iStep);
    system("python mesh_arch.py " + seed_size);


    current_directory = pwd;
    current_path = split(current_directory,DIR_DELIMINATOR);
    previous_directory = join(current_path(1:(end-1)),DIR_DELIMINATOR);
    geometry_path = previous_directory{1} + DIR_DELIMINATOR + "geometry" + DIR_DELIMINATOR + SYSTEM_NAME;

    if isfolder(geometry_path)
        rmdir(geometry_path,"s")
    end
    mkdir(geometry_path)

    data_path = previous_directory{1} + DIR_DELIMINATOR + "data";
    if isfolder(data_path)
        rmdir(data_path,"s")
    end
    mkdir(data_path)
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
    num_dofs(iStep) = num_dof;
    %-----------
    % plot_fe_mesh(geometry)
    % str = regexprep(string(num_dof),'(?<!\.\d*)\d{1,3}(?=(\d{3})+\>)','$&,');
    % sprintf(str + " dofs")

    %-----------
    %get eigenvalues
    cd(previous_directory{1})
    set_logging_level(3);
    set_visualisation_level(1);
    Model = Dynamic_System(SYSTEM_NAME,0,1:num_evals);
    cd(current_directory)

    evals = Model.reduced_eigenvalues;
    natural_frequency(:,iStep) = sqrt(evals);
end



figure;
semilogx(num_dofs,natural_frequency)