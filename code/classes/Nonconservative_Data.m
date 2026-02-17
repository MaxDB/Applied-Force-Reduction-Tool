classdef Nonconservative_Data
    properties
        damping_type
        damping_coefficients
        get_damping_matrix
        
        num_applied_forces
        force_type
        
        force_shape
        max_amplitude

        orth_force_shape
        orth_max_amplitude

        amplitude
        frequency
        phase

        Model
    end

    methods
        function obj = Nonconservative_Data(Model,Forcing_Data)
            obj.num_applied_forces = 1;
            obj.Model = Model;
            % obj = process_damping(obj,Damping_Data);
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
            max_amp = Forcing_Data.max_amplitude;
            switch obj.force_type
                case "point"
                    shape = zeros(obj.Model.num_dof,1);
                    if ~isfield(Forcing_Data,"dof_type")
                        dof_type = "pre_bc";
                    else
                        dof_type = Forcing_Data.dof_type;
                    end

                    switch dof_type
                        case "pre_bc"
                            num_dof = obj.Model.num_dof + length(obj.Model.dof_boundary_conditions);
                            shape = zeros(num_dof,1);
                            shape(Forcing_Data.dof) = 1;
                            shape(obj.Model.dof_boundary_conditions) = [];
                        case "post_bc"
                            shape(Forcing_Data.dof) = 1;
                        otherwise
                            error("")
                    end
                    
                case "shape"
                    shape = Forcing_Data.shape;
                    if size(shape,1) ~= obj.Model.num_dof
                        error("Force shape doesn't match system size")
                    end
                case "uniform"
                    num_dof = obj.Model.num_dof + + length(obj.Model.dof_boundary_conditions);
                    num_dimensions = get_num_node_dimensions(obj.Model);
                    num_nodes = num_dof/num_dimensions;
                    shape = zeros(num_dof,1);
                    direction_index = (0:(num_nodes-1))*num_dimensions + Forcing_Data.direction;
                    shape(direction_index) = 1;
                    shape(obj.Model.dof_boundary_conditions) = [];
                    shape = shape/norm(shape);
                    
                otherwise
                    Error("Unsupported force type: '" + obj.force_type + "'")
            end

            [obj.orth_force_shape,obj.orth_max_amplitude] = obj.orthogonalise_force(shape,max_amp);
            obj.force_shape = shape;
            obj.max_amplitude = max_amp;
        end
        %--------------------------------------------
        function [orth_force_shape,orth_max_amp] = orthogonalise_force(obj,force_shape,max_amp)
            evec_r = obj.Model.reduced_eigenvectors;
            if class(evec_r) == "Large_Matrix_Pointer"
                evec_r = evec_r.load();
            end
            num_r_modes = size(evec_r,2);
            phi_a = force_shape*max_amp;
            
            for iMode = 1:num_r_modes
                evec_i = evec_r(:,iMode);
                phi_a_r = evec_i*(evec_i'*obj.Model.mass*phi_a);
                phi_a = phi_a - phi_a_r;
                
            end
            if norm(phi_a) == 0
                error("No need to add force")
            end
            phi_a_amp = norm(phi_a);
            norm_test = phi_a'*obj.Model.mass*phi_a;
    
            
            orth_force_shape = phi_a/sqrt(norm_test);

            lambda_p = orth_force_shape'*force_shape*max_amp;
            % lambda_r = evec_r'*force_shape*max_amp;
            orth_max_amp = lambda_p;
        end
    end

end