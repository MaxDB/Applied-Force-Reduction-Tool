function delete_cache(system_name,type,varargin)
if isa(system_name,"Dynamic_System")
    system_name = system_name.system_name;
end

geometry_path = "geometry\" + system_name+ "\";

switch type
    case "matrices"
        file_path = geometry_path + "matrices.mat";
        if isfile(file_path)
            delete(file_path)
        end
    case "force"
        file_path = geometry_path + "force_calibration.mat";
        if ~isfile(file_path)
            return
        end
        load(file_path,"Force_Calibration")
        if nargin == 2
            delete(file_path)
            return
        end
        energy_limit = varargin{1};
        energy_index = find(Force_Calibration.energy_limit == energy_limit);

        if nargin == 3
            Force_Calibration.energy_limit(energy_index) = [];
            Force_Calibration.force_limit(energy_index) = [];
            Force_Calibration.calibrated_modes(energy_index) = [];
        else

            delete_modes = varargin{2}';
            num_energy_calibrations = length(energy_index);
            for iCalibration = 1:num_energy_calibrations
                index = energy_index(iCalibration);
                calibrated_modes = Force_Calibration.calibrated_modes{index};
                calibrated_forces = Force_Calibration.force_limit{index};
                num_modes = size(calibrated_modes,1);
                for iMode = 1:num_modes
                    jMode = (num_modes+1) - iMode;
                    if isequal(calibrated_modes(jMode),delete_modes)
                        calibrated_modes(jMode) = [];
                        calibrated_forces(jMode,:) = [];
                    end
                end

                if isempty(calibrated_modes)
                    Force_Calibration.energy_limit(index) = [];
                    Force_Calibration.force_limit(index) = [];
                    Force_Calibration.calibrated_modes(index) = []; 
                else
                    Force_Calibration.force_limit{index} = calibrated_forces;
                    Force_Calibration.calibrated_modes{index} = calibrated_modes;
                end

            end
        end
        save(file_path,"Force_Calibration")
    case "mesh_data"
        file_path = geometry_path + "mesh_data.mat";
        if isfile(file_path)
            delete(file_path)
        end
end