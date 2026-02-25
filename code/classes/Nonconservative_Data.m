classdef Nonconservative_Data
    properties
        damping_type
        damping_coefficients
        get_damping_matrix
        
        num_applied_forces
        force_type
        
        force_shape
        orth_force_shape

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

            obj.orth_force_shape = obj.orthogonalise_force(shape);
            obj.force_shape = shape;
        end
        %--------------------------------------------
        function orth_force_shape = orthogonalise_force(obj,force_shape)

            evec_r = obj.Model.reduced_eigenvectors;
            % small_disp = obj.Model.stiffness\force_shape;
            % phi_a = small_disp - evec_r*(evec_r'*obj.Model.mass*small_disp);
            applied_force = force_shape - obj.Model.mass*evec_r*(evec_r'*force_shape);
            phi_a = obj.Model.mass\applied_force;

            if norm(phi_a) == 0
                error("No need to add force")
            end
            norm_test = phi_a'*obj.Model.mass*phi_a;
            orth_force_shape = phi_a/sqrt(norm_test);
        end
    end

end