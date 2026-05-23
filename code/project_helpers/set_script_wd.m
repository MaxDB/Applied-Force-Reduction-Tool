function file_path = set_script_wd
file_path = matlab.desktop.editor.getActiveFilename;
target_wd = fileparts(file_path);
if ~isequal(target_wd,pwd)
    cd(target_wd)
    log_message = sprintf("Working directory changed to script folder");
    logger(log_message,1)
end
end