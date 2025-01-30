function reset_geometry_directory(system_name)
file_names = {system_name + ".inp","mesh_data.mat","force_calibration.mat","matrices.mat"};
geometry_path = "geometry\" + system_name + "\";

if ~isfolder(geometry_path)
    warning("'" + system_name + "' not found locally")
end

num_files = length(file_names);
for iFile = 1:num_files
    file_path = geometry_path + file_names{iFile};
    if isfile(file_path)
        delete(file_path)
    end
end
end