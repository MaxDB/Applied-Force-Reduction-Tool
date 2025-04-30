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
                if isequal(Force_Calibration.calibrated_modes{index},delete_modes)
                    Force_Calibration.energy_limit(index) = [];
                    Force_Calibration.force_limit(index) = [];
                    Force_Calibration.calibrated_modes(index) = [];
                end
            end
        end
        save(file_path,"Force_Calibration")
end