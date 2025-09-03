function system_name = get_system_name(varargin)

switch nargin
    case 1
        obj = varargin{1};
        switch class(obj)
            case "Static_Dataset"
                Model = obj.Model;
            case "Dynamic_System"
                Model = obj;
        end
    case 2
        Model.system_name = varargin{1};
        Model.reduced_modes = varargin{2};

end

r_modes = Model.reduced_modes;
mode_id = join(string(r_modes),"");
system_name = Model.system_name + "_" + mode_id; 
end