classdef Reduced_System
    %Stores the different polynomial models required for ICE-IC
    properties
        Force_Polynomial
        Physical_Displacement_Polynomial
            
        Potential_Polynomial
        Reduced_Stiffness_Polynomial

        Low_Frequency_Stiffness_Polynomial
        Low_Frequency_Coupling_Gradient_Polynomial
        
        minimum_displacement
        Model
        Dynamic_Validation_Data

        data_path
        reduced_displacement_limits

    end
    methods
        function obj = Reduced_System(Static_Data,degree)
            SHIFT_ON = 1;
            SCALE_ON = 1;
            % obj.minimum_displacement = 1e-16; %any condensed displacement with smaller range will be treated as 0
            obj.minimum_displacement = Static_Data.Model.Static_Options.minimum_displacement;

            if nargin == 1
                degree = Static_Data.verified_degree;
            end

            if isscalar(degree)
                force_degree = degree;
                disp_degree = degree;
            else
                force_degree = degree(1);
                disp_degree = degree(2);
            end
            
            obj.Model = Static_Data.Model; 

            r = Static_Data.get_dataset_values("reduced_displacement");
            f = Static_Data.get_dataset_values("restoring_force");
            displacement = Static_Data.get_dataset_values("physical_displacement");
            eval_r = Static_Data.Model.reduced_eigenvalues;
            evec_r = Static_Data.Model.reduced_eigenvectors;
            
            obj.reduced_displacement_limits = [min(r,[],2),max(r,[],2)];
           
            Force_Poly = Polynomial(r,f,force_degree,"constraint",{"linear_force",eval_r},"coupling","force","shift",SHIFT_ON,"scale",SCALE_ON);
            % Condensed_Poly = Polynomial(r,theta,disp_degree,"constraint",{"linear_disp",0},"shift",SHIFT_ON,"scale",SCALE_ON,"minimum_output",obj.MINIMUM_DISPLACEMENT);
            Displacement_Poly = Polynomial(r,displacement,disp_degree,"constraint",{"linear_disp",evec_r},"shift",SHIFT_ON,"scale",SCALE_ON,"minimum_output",obj.minimum_displacement);

            Potential_Poly = integrate_polynomial(Force_Poly);
            Reduced_Stiffness_Poly = differentiate_polynomial(Force_Poly);

            

            obj.Force_Polynomial = Force_Poly;
            obj.Physical_Displacement_Polynomial = Displacement_Poly;
            obj.Potential_Polynomial = Potential_Poly;
            obj.Reduced_Stiffness_Polynomial = Reduced_Stiffness_Poly;
            
            static_data_path = split(Static_Data.get_data_path,"\");
            obj.data_path = join(static_data_path(1:2),"\") + "\";



            if isempty(Static_Data.Dynamic_Validation_Data)
                return
            end
            
            obj.Dynamic_Validation_Data = Static_Data.Dynamic_Validation_Data;
            
            h_stiffness = Static_Data.get_dataset_values("low_frequency_stiffness");
            h_coupling_gradient = Static_Data.get_dataset_values("low_frequency_coupling_gradient");

            h_stiffness_degree = degree(3);
            h_coupling_gradient_degree = degree(4);

            h_stiffness_0 = Static_Data.Dynamic_Validation_Data.h_stiffness_0;
            h_coupling_gradient_0 = Static_Data.Dynamic_Validation_Data.h_coupling_gradient_0;

            H_Stiffness_Poly = Polynomial(r,h_stiffness,h_stiffness_degree,"constraint",{"constant",h_stiffness_0},"coupling","stiffness","shift",SHIFT_ON,"scale",SCALE_ON);
            H_Coupling_Gradient_Poly = Polynomial(r,h_coupling_gradient,h_coupling_gradient_degree,"constraint",{"constant",h_coupling_gradient_0},"shift",SHIFT_ON,"scale",SCALE_ON);
            
            obj.Low_Frequency_Stiffness_Polynomial = H_Stiffness_Poly;
            obj.Low_Frequency_Coupling_Gradient_Polynomial = H_Coupling_Gradient_Poly;
        end
        %-----------------------------------------------------------------%


        %-----------------------------------------------------------------%
        function x = expand(obj,r,varargin)
            r_evec = obj.Model.reduced_eigenvectors;
            x_Poly = obj.Physical_Displacement_Polynomial;
            x = x_Poly.evaluate_polynomial(r);
            
            if nargin == 2
                return
            end
            
            %needs to be cleaned up
            if isstring(varargin{1,1}) && varargin{1,1} == "full"
                node_map = obj.Model.node_mapping;
                dof_bc = max(node_map(:,1));
                x_bc = zeros(dof_bc,size(x,2));
                x_bc(node_map(:,1),:) = x(node_map(:,2),:);
                x = x_bc;
                return
            end
            h = varargin{1,1};

            L_evec = obj.get_current_L_eigenvectors;
            h_evec = [r_evec,L_evec];

            Theta_Hat_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;

            num_x = size(r,2);
            num_dof = size(h_evec,1);
            h_coupling = zeros(num_dof,num_x);
            for iX = 1:num_x   
                h_coupling(:,iX) = Theta_Hat_Poly.evaluate_polynomial(r(:,iX))*h(:,iX);
            end
            x = x + h_evec*h + h_coupling;
        end
        %-----------------------------------------------------------------%
        function x_dot = expand_velocity(obj,r,r_dot,varargin)
            r_evec = obj.Model.reduced_eigenvectors;
            x_Poly = obj.Physical_Displacement_Polynomial;
            x_dr_Poly = differentiate_polynomial(x_Poly);

            x_dot = x_dr_Poly.evaluate_polynomial(r)*r_dot;
            
 

            % num_x = size(r,2);
            % num_modes = size(r,1);
            % x_dr_prod = zeros(size(x_Poly,1),num_x);
            % for iX = 1:num_x
            %     x_dr_prod(:,iX) = x_dr_Poly.evaluate_polynomial(r(:,iX))*r_dot(:,iX);
            % end
            % 
            % x_dot = r_evec*r_dot + x_dr_prod;
            % 
            % if nargin == 3
            %     return
            % end
            % h = varargin{1};
            % h_dot = varargin{2};
            % 
            % L_evec = obj.get_current_L_eigenvectors;
            % h_evec = [r_evec,L_evec];
            % 
            % G_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;
            % Gdr_Poly = G_Poly.differentiate_polynomial;
            % 
            % num_x = size(r,2);
            % num_dof = size(h_evec,1);
            % G_h_dot_prod = zeros(num_dof,num_x);
            % G_dr_r_dot_h_prod = zeros(num_dof,num_x);
            % for iX = 1:num_x   
            %     G_h_dot_prod(:,iX) = G_Poly.evaluate_polynomial(r(:,iX))*h_dot(:,iX);
            %     G_dr_r_dot_h_prod(:,iX) = tensorprod(Gdr_Poly.evaluate_polynomial(r(:,iX)),r_dot(:,iX),3,1)*h(:,iX);
            % end
            % x_dot = x_dot + h_evec*h_dot + G_h_dot_prod + G_dr_r_dot_h_prod;
        end
        %-----------------------------------------------------------------%
        function beta_bar = get_beta_bar(obj,Poly_1,Poly_2)
            if nargin == 2
                Poly_2 = Poly_1;
            end
            coeffs_1 = Poly_1.coefficients;
            coeffs_2 = Poly_2.coefficients;

            mass = obj.Model.mass;
            switch ndims(coeffs_2)
                case 2
                    beta_bar = coeffs_1*mass*coeffs_2';
                case 3
                    coeffs_2 = permute(coeffs_2,[2,3,1]);
                    pre_beta_bar = tensorprod(coeffs_1,full(mass),2,1);
                    beta_bar = tensorprod(pre_beta_bar,coeffs_2,ndims(pre_beta_bar),1);
            end
        end
        %-----------------------------------------------------------------%
        function beta_mode = get_beta_mode(obj,left_product,right_product)
            mass = obj.Model.mass;
            if ndims(left_product) >= ndims(right_product)
                switch ndims(left_product)
                    case 2
                        beta_mode = left_product*mass*right_product;
                    case 3
                        pre_beta_mode = tensorprod(left_product,full(mass),2,1);
                        beta_mode = tensorprod(pre_beta_mode,right_product,ndims(pre_beta_mode),1);
                end
            else
                switch ndims(right_product)
                    case 2
                        beta_mode = left_product*mass*right_product;
                    case 3
                        pre_beta_mode = tensorprod(full(mass),right_product,2,1);
                        beta_mode = tensorprod(left_product,pre_beta_mode,2,1);
                end
            end

        end
        %-----------------------------------------------------------------%
        function input_index = get_max_input_order(obj)
            num_modes = length(obj.Model.reduced_modes);
            degree(1) = obj.Force_Polynomial.polynomial_degree;
            degree(2) = obj.Physical_Displacement_Polynomial.polynomial_degree;
            if ~isempty(obj.Dynamic_Validation_Data)
                degree(3) = obj.Low_Frequency_Stiffness_Polynomial.polynomial_degree;
                degree(4) = obj.Low_Frequency_Coupling_Gradient_Polynomial.polynomial_degree;
            end
            
            max_degree = max(degree);
            input_index = Polynomial.get_input_index(max_degree,num_modes);
        end
        %-----------------------------------------------------------------%
        function Eom_Input = get_solver_inputs(obj,type,varargin)
            switch type
                case "coco_backbone"
                    input_order = obj.get_max_input_order;

                    Force_Data.coeffs = obj.Force_Polynomial.coefficients';
                    Force_Data.scale_factor = obj.Force_Polynomial.scaling_factor;
                    Force_Data.shift_factor = obj.Force_Polynomial.shifting_factor;
                    Force_Diff_Data = obj.Force_Polynomial.get_diff_data(1);
                    Force_Data.diff_scale_factor = Force_Diff_Data.diff_scale_factor;
                    Force_Data.diff_mapping = Force_Diff_Data.diff_mapping;

                    Disp_Data.beta_bar = obj.get_beta_bar(obj.Physical_Displacement_Polynomial);
                    Disp_Data.scale_factor = obj.Physical_Displacement_Polynomial.scaling_factor;
                    Disp_Data.shift_factor = obj.Physical_Displacement_Polynomial.shifting_factor;
                    Disp_Diff_Data = obj.Physical_Displacement_Polynomial.get_diff_data(3);
                    Disp_Data.diff_scale_factor = Disp_Diff_Data.diff_scale_factor;
                    Disp_Data.diff_mapping = Disp_Diff_Data.diff_mapping;


                    Eom_Input.input_order = input_order;
                    Eom_Input.Force_Data = Force_Data;
                    Eom_Input.Disp_Data = Disp_Data;
                    Eom_Input.Potential_Polynomial = obj.Potential_Polynomial;
                    Eom_Input.energy_limit = obj.Model.energy_limit;
                case "h_prediction"
                    input_order = obj.get_max_input_order;
                    scale_factor = obj.Force_Polynomial.scaling_factor;
                    shift_factor = obj.Force_Polynomial.shifting_factor;
                    
                    Reduced_Force_Data.coeffs = obj.Force_Polynomial.coefficients';
                     
                    Physical_Disp_Diff_Data = obj.Physical_Displacement_Polynomial.get_diff_data(2);
                    Physical_Disp_Data.diff_scale_factor = Physical_Disp_Diff_Data.diff_scale_factor;
                    Physical_Disp_Data.diff_mapping = Physical_Disp_Diff_Data.diff_mapping;
                    physical_disp_coeffs = obj.Physical_Displacement_Polynomial.coefficients';
                    
                    H_Stiffness_Poly = obj.Low_Frequency_Stiffness_Polynomial;
                    H_Force_Data.coeffs = permute(H_Stiffness_Poly.coefficients,[2,3,1]);
                    
                    H_Disp_Grad_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;
                    h_disp_coeff = H_Disp_Grad_Poly.coefficients;

                    Beta_Bar_Data.h_disp = obj.get_beta_mode(h_disp_coeff,permute(h_disp_coeff,[2,3,1]));
                    Beta_Bar_Data.h_disp_r_disp = obj.get_beta_mode(h_disp_coeff,physical_disp_coeffs);

                   
                    Eom_Input.input_order = input_order;
                    Eom_Input.scale_factor = scale_factor;
                    Eom_Input.shift_factor = shift_factor;
                    
                    Eom_Input.Reduced_Force_Data = Reduced_Force_Data;
                    Eom_Input.H_Force_Data = H_Force_Data;
                    Eom_Input.Physical_Disp_Data = Physical_Disp_Data;
                    Eom_Input.Beta_Bar_Data = Beta_Bar_Data;

                case "h_analysis"
                    input_order = obj.get_max_input_order;
                    scale_factor = obj.Force_Polynomial.scaling_factor;
                    shift_factor = obj.Force_Polynomial.shifting_factor;

                    Physical_Disp_Diff_Data = obj.Physical_Displacement_Polynomial.get_diff_data(1);
                    Physical_Disp_Data.diff_scale_factor = Physical_Disp_Diff_Data.diff_scale_factor;
                    Physical_Disp_Data.diff_mapping = Physical_Disp_Diff_Data.diff_mapping;
                    Physical_Disp_Data.Poly = obj.Physical_Displacement_Polynomial;
                    physical_disp_coeffs = obj.Physical_Displacement_Polynomial.coefficients';

                    Potential_Poly = obj.Potential_Polynomial;

                    H_Stiff_Poly = obj.Low_Frequency_Stiffness_Polynomial;
                    Force_Poly = obj.Force_Polynomial;
                    
                    H_Disp_Grad_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;
                    h_disp_coeff = H_Disp_Grad_Poly.coefficients;

                    L_evecs = obj.get_current_L_eigenvectors;
                    mass = obj.Model.mass;
                    L_disp_transform = (L_evecs'*mass);


                    Beta_Bar_Data.h_disp = obj.get_beta_mode(h_disp_coeff,permute(h_disp_coeff,[2,3,1]));
                    Beta_Bar_Data.h_disp_r_disp = obj.get_beta_mode(h_disp_coeff,physical_disp_coeffs);
                    Beta_Bar_Data.r_disp = obj.get_beta_mode(physical_disp_coeffs',physical_disp_coeffs);

                    Eom_Input.input_order = input_order;
                    Eom_Input.scale_factor = scale_factor;
                    Eom_Input.shift_factor = shift_factor;
                    
                    Eom_Input.H_Force_Poly = H_Stiff_Poly;
                    Eom_Input.Physical_Disp_Data = Physical_Disp_Data;
                    Eom_Input.Potential_Poly = Potential_Poly;
                    Eom_Input.Force_Poly = Force_Poly;
                    % Eom_Input.H_Disp_Data = H_Disp_Data;
                    Eom_Input.Beta_Bar_Data = Beta_Bar_Data;
                    Eom_Input.L_disp_transform = L_disp_transform;

                case "coco_frf"
                    Eom_Input = obj.get_solver_inputs("coco_backbone");
                    Nc_Inputs = varargin{1,1};
                    
                    displacement_coeffs = obj.Physical_Displacement_Polynomial.coefficients';
                    damping_beta = displacement_coeffs'*Nc_Inputs.damping*displacement_coeffs;
                    Eom_Input.Damping_Data.damping_beta = damping_beta;

                    switch Nc_Inputs.force_type
                        case "modal"
                            Eom_Input.Applied_Force_Data.shape = @(t,amp,T) modal_force(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dx = @(t,amp,T) modal_force_dx(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dTper = @(t,amp,T) modal_force_dTper(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dA = @(t,amp,T) modal_force_dA(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dt = @(t,amp,T) modal_force_dt(t,amp,T,Nc_Inputs.mode_map);

                            Eom_Input.Applied_Force_Data.type = Nc_Inputs.force_type;
                            switch Nc_Inputs.continuation_variable
                                case "amplitude"
                                    Eom_Input.Applied_Force_Data.period = 2*pi/Nc_Inputs.frequency;
                                case "frequency"
                                    Eom_Input.Applied_Force_Data.amplitude = Nc_Inputs.amplitude;
                            end
                        case "point force"
                            num_r_modes = size(obj.Model.reduced_modes,2);

                            Eom_Input.Applied_Force_Data.shape = @(t,amp,T) sine_force(t,amp,T);
                            Eom_Input.Applied_Force_Data.shape_dx = @(t,amp,T) sine_force_dx(t,amp,T,num_r_modes);
                            Eom_Input.Applied_Force_Data.shape_dTper = @(t,amp,T) sine_force_dTper(t,amp,T);
                            Eom_Input.Applied_Force_Data.shape_dA = @(t,amp,T) sine_force_dA(t,amp,T);
                            Eom_Input.Applied_Force_Data.shape_dt = @(t,amp,T) sine_force_dt(t,amp,T);

                            Eom_Input.Applied_Force_Data.type = Nc_Inputs.force_type;
                            switch Nc_Inputs.continuation_variable
                                case "amplitude"
                                    Eom_Input.Applied_Force_Data.period = 2*pi/Nc_Inputs.frequency;
                                case "frequency"
                                    Eom_Input.Applied_Force_Data.amplitude = Nc_Inputs.amplitude;
                            end

                            h_disp_force_beta = displacement_coeffs'*Nc_Inputs.amplitude_shape;
                            Eom_Input.Applied_Force_Data.disp_force_beta = h_disp_force_beta;

                        otherwise 
                            error("Unknown force type: '" + Nc_Inputs.force_type + "'")
                    end
                case "forced_h_prediction"
                    Eom_Input = obj.get_solver_inputs("h_prediction");
                    Nc_Inputs = varargin{1,1};

                    damping = Nc_Inputs.damping;

                    displacement_coeffs = obj.Physical_Displacement_Polynomial.coefficients';

                    h_disp_coeff = obj.Low_Frequency_Coupling_Gradient_Polynomial.coefficients;
                    h_disp_coeff_T = permute(h_disp_coeff,[2,3,1]);
                    damping_h_disp = tensorprod(full(damping),h_disp_coeff_T,2,1);
                    Beta_Damping.h_disp = obj.get_beta_mode(h_disp_coeff,damping_h_disp);
                    Beta_Damping.h_disp_r_disp = obj.get_beta_mode(h_disp_coeff,damping*displacement_coeffs);
 
                    Eom_Input.Beta_Damping = Beta_Damping;


                    switch Nc_Inputs.force_type
                        case "modal"
                            Eom_Input.Applied_Force_Data.shape = @(t,amp,T) modal_force(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dx = @(t,amp,T) modal_force_dx(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dTper = @(t,amp,T) modal_force_dTper(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dA = @(t,amp,T) modal_force_dA(t,amp,T,Nc_Inputs.mode_map);
                            Eom_Input.Applied_Force_Data.shape_dt = @(t,amp,T) modal_force_dt(t,amp,T,Nc_Inputs.mode_map);

                            Eom_Input.Applied_Force_Data.type = Nc_Inputs.force_type;
                            switch Nc_Inputs.continuation_variable
                                case "amplitude"
                                    Eom_Input.Applied_Force_Data.period = 2*pi/Nc_Inputs.frequency;
                                case "frequency"
                                    Eom_Input.Applied_Force_Data.amplitude = Nc_Inputs.amplitude;
                            end
                        case "point force"
                            num_r_modes = size(obj.Model.reduced_modes,2);

                            Eom_Input.Applied_Force_Data.shape = @(t,amp,T) sine_force(t,amp,T);
                            Eom_Input.Applied_Force_Data.shape_dx = @(t,amp,T) sine_force_dx(t,amp,T,num_r_modes);
                            Eom_Input.Applied_Force_Data.shape_dTper = @(t,amp,T) sine_force_dTper(t,amp,T);
                            Eom_Input.Applied_Force_Data.shape_dA = @(t,amp,T) sine_force_dA(t,amp,T);
                            Eom_Input.Applied_Force_Data.shape_dt = @(t,amp,T) sine_force_dt(t,amp,T);

                            Eom_Input.Applied_Force_Data.type = Nc_Inputs.force_type;
                            switch Nc_Inputs.continuation_variable
                                case "amplitude"
                                    Eom_Input.Applied_Force_Data.period = 2*pi/Nc_Inputs.frequency;
                                case "frequency"
                                    Eom_Input.Applied_Force_Data.amplitude = Nc_Inputs.amplitude;
                            end

                            h_disp_force_beta = obj.get_beta_mode(h_disp_coeff,Nc_Inputs.amplitude_shape);
                            Eom_Input.Applied_Force_Data.h_disp_force_beta = h_disp_force_beta;

                        otherwise
                            error("Unknown force type: '" + Nc_Inputs.force_type + "'")
                    end
                case "forced_h_analysis"
                    Eom_Input = obj.get_solver_inputs("h_analysis");
                    Nc_Inputs = varargin{1,1};
                   
            end
        end
        %-----------------------------------------------------------------%
        function L_eigenvectors = get_current_L_eigenvectors(obj)
            L_modes = obj.Dynamic_Validation_Data.current_L_modes;
            all_L_modes = obj.Model.low_frequency_modes;
            [~,L_map] = ismember(L_modes,all_L_modes);
            L_eigenvectors = obj.Model.low_frequency_eigenvectors(:,L_map);
        end
        %-----------------------------------------------------------------%
    end

end