function project_path = get_project_path
FILE_SEPERATOR = "\";
working_directory = pwd;
path_components = split(working_directory,FILE_SEPERATOR);
path_end = find(path_components == "Applied Force Reduction");
project_path = join(path_components(1:path_end),FILE_SEPERATOR);
project_path = project_path{1};
end