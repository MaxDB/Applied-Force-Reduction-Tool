function delete_dynamic_data(system_name)
Dyn_Data = initalise_dynamic_data(system_name);
if isempty(Dyn_Data)
    return
end

Dyn_Data.remove_solutions("all");
delete("data\" + system_name + "\Dyn_Data.mat")
end