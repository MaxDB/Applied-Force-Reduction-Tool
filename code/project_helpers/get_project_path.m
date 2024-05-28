function project_path = get_project_path
FILE_SEPERATOR = "\";
working_directory = pwd;
path_components = split(working_directory,FILE_SEPERATOR);
project_path = join(path_components(1:(end-1)),FILE_SEPERATOR);
project_path = project_path{1};
end