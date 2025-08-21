function system_name = get_system_name(obj)
switch class(obj)
    case "Static_Dataset"
        Model = obj.Model;
    case "Dynamic_System"
        Model = obj;
end

r_modes = Model.reduced_modes;
mode_id = join(string(r_modes),"");
system_name = Model.system_name + "_" + mode_id; 
end