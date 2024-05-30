function Dyn_Data = initalise_dynamic_data(system_name,Additional_Output)
file_path = "data\" + system_name + "\";
if isfile(file_path + "Dyn_Data.mat")
    load(file_path + "Dyn_Data.mat","Dyn_Data");
else
    load(file_path + "Static_Data.mat","Static_Data");
    Rom = Reduced_System(Static_Data);
    if nargin == 1
        Additional_Output.type = "none";
    end
    Dyn_Data = Dynamic_Data(Rom,Additional_Output);
end
end