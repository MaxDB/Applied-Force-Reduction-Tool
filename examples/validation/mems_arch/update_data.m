clear
close all
num_workers = 6;

old_index = 1;
new_index = 1;

data_path = "data\size_data\workers_" + num_workers + "\";

old_data = "size_data";
new_data = "size_data_new";

size_data_old = load(data_path + old_data,"Size_Data");
size_data_new = load(data_path + new_data,"Size_Data");

Size_Data_Old = size_data_old.Size_Data;
Size_Data_New = size_data_new.Size_Data;

if isfield(Size_Data_Old,"dynamic_time")
    Size_Data_Old = rmfield(Size_Data_Old,"dynamic_time");
end

Size_Data_Old(old_index) = Size_Data_New(new_index);
Size_Data = Size_Data_Old;
save(data_path + old_data,"Size_Data")