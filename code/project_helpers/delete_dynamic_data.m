function delete_dynamic_data(system_name)
file_name = "data\" + system_name + "\Dyn_Data.mat";
if ~isfile(file_name)
    return
end


Dyn_Data = initalise_dynamic_data(system_name);
if isempty(Dyn_Data)
    return
end

Dyn_Data.remove_solution("all")
delete(file_name)

end