function Dyn_Data = initalise_dynamic_data(system_name)
file_path = "data\" + system_name + "\static_data\";
if isfile(file_path + "Dyn_Data.mat")
    load(file_path + "Dyn_Data.mat","Dyn_Data");
else
    load(file_path + "Static_Data.mat","Static_Data");
    Rom = Reduced_System(Static_Data);
    Dyn_Data = Dynamic_Data(Rom);
end
end