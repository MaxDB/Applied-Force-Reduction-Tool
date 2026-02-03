function Static_Data = load_static_data(system_name,varargin)
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);


load_nc = 0;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "nonconservative"
            load_nc = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------%

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
    file_path = "data\" + system_name + "\static_data";
    if load_nc
        file_path = file_path + "_nc";
    end
    file_path = file_path + "\Static_Data.mat";
    
    if isfile(file_path)
        load(file_path,"Static_Data")
    else
        Static_Data = [];
    end
end