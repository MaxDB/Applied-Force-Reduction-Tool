function validation_modes = process_validation_modes(mode_list,Model)

if isstring(mode_list)
    switch mode_list
        case "all"
            switch Model.system_type
                case "indirect"
                    validation_modes = Model.low_frequency_modes;
                case "direct"
                    validation_modes = 1:Model.num_dof;
            end
        otherwise
            num_modes = size(mode_list,2);
            validation_modes = zeros(1,num_modes);
            for iMode = 1:num_modes
                mode = mode_list(iMode);
                mode_num = extract(mode,digitsPattern);
                mode_type = extract(mode,lettersPattern);
                validation_modes(1,iMode) = double(mode_num);
                if mode_type == "f"
                    validation_modes(1,iMode) = validation_modes(1,iMode) + 2000;
                end
            end
    end
else
    validation_modes = mode_list;
end

end