function project_path = get_project_path
FILE_SEPERATOR = "\";
function_path = mfilename("fullpath");

path_components = split(function_path,FILE_SEPERATOR);
project_path = join(path_components(1:(end-3)),FILE_SEPERATOR);
project_path = project_path{1};
end