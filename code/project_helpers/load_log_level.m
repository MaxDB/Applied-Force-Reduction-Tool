function log_level = load_log_level(type)
DEFAULT_LEVEL = 1;

switch type
    case "log"
        var_name = "logging_level";
        log_path = "data\log_level.mat";
    case "plot"
        var_name = "plotting_level";
        log_path = "data\plot_level.mat";
end

if ~isfile(log_path)
    log_level = DEFAULT_LEVEL;
    return
end

var_struct = load(log_path,var_name);
log_level = var_struct.(var_name);

end