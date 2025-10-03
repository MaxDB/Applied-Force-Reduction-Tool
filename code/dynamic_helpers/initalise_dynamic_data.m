function Dyn_Data = initalise_dynamic_data(system_name)
initalisation_time_start = tic;
Dyn_Data = init_helper(system_name);
initalisation_time = toc(initalisation_time_start);
log_message = sprintf("Dynamic Dataset initalised: %.1f seconds" ,initalisation_time);
logger(log_message,1)
end

function Dyn_Data = init_helper(system_name)
if class(system_name) == "Dynamic_Dataset"
    Dyn_Data = system_name;
    return
end
file_path = "data\" + system_name + "\";
if isfile(file_path + "Dyn_Data.mat")
    load(file_path + "Dyn_Data.mat","Dyn_Data");
else
    Static_Data = load_static_data(system_name);
    if isempty(Static_Data)
        Dyn_Data = [];
        return
    end
    Rom = Reduced_System(Static_Data);
    Dyn_Data = Dynamic_Dataset(Rom);
end
end