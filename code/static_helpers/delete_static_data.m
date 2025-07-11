function delete_static_data(system_name)
switch class(system_name)
    case "Dynamic_Dataset"
        Model = system_name.Dynamic_Model.Model;
        r_modes = Model.reduced_modes;
        mode_id = join(string(r_modes),"");

        system_name = Model.system_name + "_" + mode_id;
    case "Reduced_System"
        Model = system_name.Model;
        r_modes = Model.reduced_modes;
        mode_id = join(string(r_modes),"");

        system_name = Model.system_name + "_" + mode_id;
    case "string"

end
    % file_path = "data\" + system_name + "\static_data";
    if endsWith(system_name,"_" + digitsPattern)
        file_path = "data\" + system_name;
    else
        file_path = "data";
    end
    if isfolder(file_path)
        rmdir(file_path,"s")
    end
end