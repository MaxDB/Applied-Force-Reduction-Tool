classdef Nonconservative_Data
    properties
        damping_type
        damping_coefficients
        get_damping_matrix
        
        num_applied_forces
        force_type
        force_shape

        max_amplitude
        frequency_range

        amplitude
        frequency
        phase

        Model
    end

    methods
        function obj = Nonconservative_Data(Model,Damping_Data,Forcing_Data)
            obj.num_applied_forces = 1;
            obj.Model = Model;
            obj = process_damping(obj,Damping_Data);
            obj = process_forcing(obj,Forcing_Data);

        end
        %--------------------------------------------
        function obj = process_damping(obj,Damping_Data)
            obj.damping_type = lower(Damping_Data.type);
            switch obj.damping_type
                case "rayleigh"
                    coeffs = [Damping_Data.mass_factor,Damping_Data.stiffness_factor];
                    obj.damping_coefficients = coeffs;
                    obj.get_damping_matrix = @() Damping_Data.mass_factor*obj.Model.mass + Damping_Data.stiffness_factor*obj.Model.stiffness;
                otherwise
                    Error("Unsupported damping type: '" + obj.damping_type + "'")
                
            end

        end
        %--------------------------------------------
        function obj = process_forcing(obj,Forcing_Data)
            obj.force_type = Forcing_Data.type;
            obj.max_amplitude = Forcing_Data.max_amplitude;
            switch obj.force_type
                case "point"
                    shape = zeros(obj.Model.num_dof,1);
                    shape(Forcing_Data.dof) = 1;
                    obj.force_shape = shape;
                otherwise
                    Error("Unsupported force type: '" + obj.force_type + "'")
            end

        end
    end

end