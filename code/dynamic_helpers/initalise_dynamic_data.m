function Dyn_Data = initalise_dynamic_data(system_name)
if class(system_name) == "Dynamic_Dataset"
    Dyn_Data = system_name;
    return
end
file_path = "data\" + system_name + "\";
if isfile(file_path + "Dyn_Data.mat")
    load(file_path + "Dyn_Data.mat","Dyn_Data");
else
    Static_Data = load_static_data(system_name);
    Rom = Reduced_System(Static_Data);
    Dyn_Data = Dynamic_Dataset(Rom);
end
end