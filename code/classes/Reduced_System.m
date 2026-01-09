classdef Reduced_System
    %Stores the different polynomial models required for the applied force
    %method
    properties
        Force_Polynomial
        Physical_Displacement_Polynomial
            
        Potential_Polynomial
        Reduced_Stiffness_Polynomial

        Low_Frequency_Stiffness_Polynomial
        Low_Frequency_Coupling_Gradient_Polynomial
        get_low_frequency_displacement
        
        minimum_displacement
        Model
        Dynamic_Validation_Data

        data_path
        reduced_displacement_limits

        id

    end
    methods
        function obj = Reduced_System(Static_Data,varargin)
            SHIFT_ON = 1;
            SCALE_ON = 1;
            obj.minimum_displacement = Static_Data.Model.Static_Options.minimum_displacement;
            
            Experimental_Opts = read_default_options("experimental"); 
            param_style = Experimental_Opts.param_style;

            %Optional argumanents
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            degree = Static_Data.verified_degree;
            rom_id = 1;
            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "degree"
                        degree = keyword_values{arg_counter};
                    case "id"
                        rom_id = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %--------------------------------
    


            if isscalar(degree)
                force_degree = degree;
                disp_degree = degree;
            else
                force_degree = degree(1);
                disp_degree = degree(2);
            end
            %-----------------




            obj.Model = Static_Data.Model; 
            obj.id = rom_id;
            
           
            f = Static_Data.get_dataset_values("restoring_force");
            displacement = Static_Data.get_dataset_values("physical_displacement");
            eval_r = Static_Data.Model.reduced_eigenvalues;
            evec_r = Static_Data.Model.reduced_eigenvectors;
            evec_r = load_data(evec_r);
            
 
           

            switch param_style
                case "disp"
                    r = Static_Data.get_dataset_values("reduced_displacement");
                    manifold_param = r;
                case "force"
                    manifold_param = f./eval_r;
            end
            obj.reduced_displacement_limits = [min(manifold_param,[],2),max(manifold_param,[],2)];

            Force_Poly = Polynomial(manifold_param,f,force_degree,"constraint",{"linear_force",eval_r},"coupling","force","shift",SHIFT_ON,"scale",SCALE_ON);
            Displacement_Poly = Polynomial(manifold_param,displacement,disp_degree,"constraint",{"linear_disp",evec_r},"shift",SHIFT_ON,"scale",SCALE_ON,"minimum_output",obj.minimum_displacement);
     
            
            % ax = plot_static_data("force",Static_Data);
            % ax = Force_Poly.plot_polynomial("axes",ax);
            % % 
            % ax = plot_static_data("force",Static_Data);
            % ax = Force_Poly.plot_polynomial("axes",ax);
            
            switch param_style
                case "disp"
                    Potential_Poly = integrate_polynomial(Force_Poly);
                    Reduced_Stiffness_Poly = differentiate_polynomial(Force_Poly);
                case "force"
                    V = Static_Data.get_dataset_values("potential_energy");
                    Potential_Poly = Polynomial(manifold_param,V,force_degree,"constraint",{"constant",0},"shift",SHIFT_ON,"scale",SCALE_ON);
                    Reduced_Stiffness_Poly = [];
            end

            % ax = plot_static_data("potential",Static_Data);
            % ax = Potential_Poly.plot_polynomial("axes",ax);
            % Reduced_Stiffness_Poly.plot_polynomial("outputs",[1,1]);
            
            

            obj.Force_Polynomial = Force_Poly;
            obj.Physical_Displacement_Polynomial = Displacement_Poly;
            obj.Potential_Polynomial = Potential_Poly;
            obj.Reduced_Stiffness_Polynomial = Reduced_Stiffness_Poly;
            
            static_data_path = split(Static_Data.get_data_path,"\");
            obj.data_path = join(static_data_path(1:2),"\") + "\";

            rom_data_path = obj.data_path + "rom_data_" + obj.id;
            if isfolder(rom_data_path)
                rmdir(rom_data_path,"s")
            end
            mkdir(rom_data_path)

            if isempty(Static_Data.Dynamic_Validation_Data)
                return
            end
            
            obj.Dynamic_Validation_Data = Static_Data.Dynamic_Validation_Data;
            
            h_stiffness = Static_Data.get_dataset_values("low_frequency_stiffness");
            h_coupling_gradient = Static_Data.get_dataset_values("low_frequency_coupling_gradient");

            h_stiffness_degree = Static_Data.Dynamic_Validation_Data.degree(1);
            h_coupling_gradient_degree = Static_Data.Dynamic_Validation_Data.degree(2);

            h_stiffness_0 = Static_Data.Dynamic_Validation_Data.h_stiffness_0;
            h_coupling_gradient_0 = Static_Data.Dynamic_Validation_Data.h_coupling_gradient_0;

            H_Stiffness_Poly = Polynomial(r,h_stiffness,h_stiffness_degree,"constraint",{"constant",h_stiffness_0},"coupling","stiffness","shift",SHIFT_ON,"scale",SCALE_ON);
            H_Coupling_Gradient_Poly = Polynomial(r,h_coupling_gradient,h_coupling_gradient_degree,"constraint",{"constant",h_coupling_gradient_0},"shift",SHIFT_ON,"scale",SCALE_ON);
            
            obj.Low_Frequency_Stiffness_Polynomial = H_Stiffness_Poly;
            obj.Low_Frequency_Coupling_Gradient_Polynomial = H_Coupling_Gradient_Poly;

            obj.get_low_frequency_displacement = obj.create_modal_validation_polynomial;
        end
        %-----------------------------------------------------------------%


        %-----------------------------------------------------------------%
        function x = expand(obj,r,varargin)
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

            H_Grad_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;

            num_x = size(r,2);
            num_dof = size(x,1);
            h_coupling = zeros(num_dof,num_x);
            for iX = 1:num_x   
                h_coupling(:,iX) = H_Grad_Poly.evaluate_polynomial(r(:,iX))*h(:,iX);
            end
            x = x + h_coupling;
        end
        %-----------------------------------------------------------------%
        function x_dot = expand_velocity(obj,r,r_dot,varargin)
            x_Poly = obj.Physical_Displacement_Polynomial;
            x_dr_Poly = differentiate_polynomial(x_Poly);
            
            num_time_points = size(r,2);
            num_dof = obj.Model.num_dof;
            x_dot = zeros(num_dof,num_time_points);
            for iT = 1:num_time_points
                x_dot(:,iT) = x_dr_Poly.evaluate_polynomial(r(:,iT))*r_dot(:,iT);
            end


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
        function x_dot = get_physical_velocity(obj,r,r_dot,output_dofs)
            Phy_Poly = obj.Physical_Displacement_Polynomial;
            Phy_Poly = Phy_Poly.subpoly(output_dofs);
            Phy_dr_Poly = differentiate_polynomial(Phy_Poly);
            
            
            num_points = size(r,2);
            num_outputs = size(output_dofs,2);
            x_dot = zeros(num_outputs,num_points);
            for iPoint = 1:num_points
                r_i = r(:,iPoint);
                r_dot_i = r_dot(:,iPoint);
                x_dr_i = Phy_dr_Poly.evaluate_polynomial(r_i);
                x_dot(:,iPoint) = x_dr_i*r_dot_i; 
            end
            
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
                    beta_bar = coeffs_1'*mass*coeffs_2;
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
        function Validation_Beta_Bar_Data = get_h_beta_bar(obj,h_disp_coeff_diff,physical_disp_coeffs)

            %reformat h_disp_coeffs to input differential form
            h_coeff_size = size(h_disp_coeff_diff);
            num_coeffs = h_coeff_size(3);
            num_dof = h_coeff_size(1);
            num_validation_modes = h_coeff_size(2);
            
            num_combined_coeffs = num_coeffs*num_validation_modes;
            h_disp_coeff = zeros(num_dof,num_combined_coeffs);
            for iMode = 1:num_validation_modes
                mode_index = iMode:num_validation_modes:num_combined_coeffs;
                h_disp_coeff(:,mode_index) = h_disp_coeff_diff(:,iMode,:);
            end

            h_disp_h_disp = obj.get_beta_mode(h_disp_coeff',h_disp_coeff);
            h_disp_r_disp = obj.get_beta_mode(h_disp_coeff',physical_disp_coeffs);
            
            num_physical_disp_coeffs = size(physical_disp_coeffs,2);
            h_disp_h_disp_diff = zeros(num_coeffs,num_validation_modes,num_validation_modes,num_coeffs);
            h_disp_r_disp_diff = zeros(num_coeffs,num_validation_modes,num_physical_disp_coeffs);
            for iMode = 1:num_validation_modes
                mode_index_one = iMode:num_validation_modes:num_combined_coeffs;
                h_disp_r_disp_diff(:,iMode,:) = h_disp_r_disp(mode_index_one,:);
                for jMode = 1:num_validation_modes
                    mode_index_two = jMode:num_validation_modes:num_combined_coeffs;
                    h_disp_h_disp_diff(:,iMode,jMode,:) = h_disp_h_disp(mode_index_one,mode_index_two);
                end
            end
            Validation_Beta_Bar_Data.h_disp = h_disp_h_disp_diff;
            Validation_Beta_Bar_Data.h_disp_r_disp = h_disp_r_disp_diff;
            % Beta_Bar_Data.h_disp = obj.get_beta_mode(h_disp_coeff_diff,permute(h_disp_coeff_diff,[2,3,1]));
            % Beta_Bar_Data.h_disp_r_disp = obj.get_beta_mode(h_disp_coeff_diff,physical_disp_coeffs);
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
            pre_dynamic_time_start = tic;

            %-------------------------------------------------------------------------%
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            rom_type = "standard";

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "rom_type"
                        rom_type = keyword_values{arg_counter};
                    case "additional_output"
                        Additional_Output = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                        %will need to fix later
                end
            end
            %-------------------------------------------------------------------------%


            file_name = "rom_data_" + obj.id + "\" + type;
            if rom_type ~= "standard"
                file_name = file_name + "_rom_type";
            end
            file_name = file_name + ".mat";
            rom_data = obj.data_path + file_name;
   
            if isfile(rom_data)
                load(rom_data,"Eom_Input");
                pre_dynamic_time = toc(pre_dynamic_time_start);
                log_message = sprintf("EoM loaded: %.1f seconds" ,pre_dynamic_time);
                logger(log_message,3)
                return
            end

            switch type
                case "coco_backbone"
                    input_order = obj.get_max_input_order;

                    Force_Data.coeffs = obj.Force_Polynomial.coefficients;
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

                    Reduced_Force_Data.coeffs = obj.Force_Polynomial.coefficients;

                    Physical_Disp_Diff_Data = obj.Physical_Displacement_Polynomial.get_diff_data(2);
                    Physical_Disp_Data.diff_scale_factor = Physical_Disp_Diff_Data.diff_scale_factor;
                    Physical_Disp_Data.diff_mapping = Physical_Disp_Diff_Data.diff_mapping;
                    physical_disp_coeffs = obj.Physical_Displacement_Polynomial.coefficients;

                    H_Stiffness_Poly = obj.Low_Frequency_Stiffness_Polynomial;
                    H_Force_Data.coeffs = H_Stiffness_Poly.coefficients;

                    H_Disp_Grad_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;
                    h_disp_coeff = H_Disp_Grad_Poly.coefficients;

                    Disp_Grad_Diff_Data = H_Stiffness_Poly.get_diff_data(2);
                    Disp_Grad_Data.diff_scale_factor = Disp_Grad_Diff_Data.diff_scale_factor;
                    Disp_Grad_Data.diff_mapping = Disp_Grad_Diff_Data.diff_mapping;

                    Beta_Bar_Data = obj.get_h_beta_bar(h_disp_coeff,physical_disp_coeffs);

                    % Beta_Bar_Data.h_disp = obj.get_beta_mode(h_disp_coeff,permute(h_disp_coeff,[2,3,1]));
                    % Beta_Bar_Data.h_disp_r_disp = obj.get_beta_mode(h_disp_coeff,physical_disp_coeffs);


                    Eom_Input.input_order = input_order;
                    Eom_Input.scale_factor = scale_factor;
                    Eom_Input.shift_factor = shift_factor;

                    Eom_Input.Reduced_Force_Data = Reduced_Force_Data;
                    Eom_Input.H_Force_Data = H_Force_Data;
                    Eom_Input.Physical_Disp_Data = Physical_Disp_Data;
                    Eom_Input.Disp_Grad_Data = Disp_Grad_Data;
                    Eom_Input.Beta_Bar_Data = Beta_Bar_Data;

                case "h_analysis"
                    input_order = obj.get_max_input_order;
                    scale_factor = obj.Force_Polynomial.scaling_factor;
                    shift_factor = obj.Force_Polynomial.shifting_factor;

                    Physical_Disp_Diff_Data = obj.Physical_Displacement_Polynomial.get_diff_data(1);
                    Physical_Disp_Data.diff_scale_factor = Physical_Disp_Diff_Data.diff_scale_factor;
                    Physical_Disp_Data.diff_mapping = Physical_Disp_Diff_Data.diff_mapping;
                    physical_disp_coeffs = obj.Physical_Displacement_Polynomial.coefficients;

                    Potential_Poly = obj.Potential_Polynomial;

                    H_Stiff_Poly = obj.Low_Frequency_Stiffness_Polynomial;
                    Force_Poly = obj.Force_Polynomial;

                    H_Disp_Grad_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;
                    h_disp_coeff = H_Disp_Grad_Poly.coefficients;

                    L_evecs = obj.get_current_L_eigenvectors;
                    mass = obj.Model.mass;
                    L_disp_transform = (L_evecs'*mass);

                    Beta_Bar_Data = obj.get_h_beta_bar(h_disp_coeff,physical_disp_coeffs);
                    % Beta_Bar_Data.h_disp = obj.get_beta_mode(h_disp_coeff,permute(h_disp_coeff,[2,3,1]));
                    % Beta_Bar_Data.h_disp_r_disp = obj.get_beta_mode(h_disp_coeff,physical_disp_coeffs);
                    Beta_Bar_Data.r_disp = obj.get_beta_mode(physical_disp_coeffs',physical_disp_coeffs);
                    
                    H_Stiffness_Poly = obj.Low_Frequency_Stiffness_Polynomial;
                    Disp_Grad_Diff_Data = H_Stiffness_Poly.get_diff_data(2);
                    Disp_Grad_Data.diff_scale_factor = Disp_Grad_Diff_Data.diff_scale_factor;
                    Disp_Grad_Data.diff_mapping = Disp_Grad_Diff_Data.diff_mapping;

                    Eom_Input.input_order = input_order;
                    Eom_Input.scale_factor = scale_factor;
                    Eom_Input.shift_factor = shift_factor;

                    Eom_Input.H_Force_Poly = H_Stiff_Poly;
                    Eom_Input.Physical_Disp_Data = Physical_Disp_Data;
                    Eom_Input.Potential_Poly = Potential_Poly;
                    Eom_Input.Force_Poly = Force_Poly;
                    % Eom_Input.H_Disp_Data = H_Disp_Data;
                    Eom_Input.Beta_Bar_Data = Beta_Bar_Data;
                    % Eom_Input.L_disp_transform = L_disp_transform;
                    Eom_Input.lf_disp_func = obj.get_low_frequency_displacement;
                    Eom_Input.Disp_Grad_Data = Disp_Grad_Data;

           
                    Eom_Input.Additional_Output = Additional_Output;
                    switch Additional_Output.output
                        case "physical displacement"
                            G_Disp_Poly = obj.Low_Frequency_Coupling_Gradient_Polynomial;
                            Eom_Input.Add_Output_Data.G_Grad_Subpoly = G_Disp_Poly.subpoly(Additional_Output.control_dof);
                            Eom_Input.Add_Output_Data.x_Subpoly = obj.Physical_Displacement_Polynomial.subpoly(Additional_Output.control_dof);
                    end

                case "coco_frf"
                    Eom_Input = obj.get_solver_inputs("coco_backbone");
                    Nc_Inputs = varargin{1,1};

                    displacement_coeffs = obj.Physical_Displacement_Polynomial.coefficients;
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
                        case "none"
                            Eom_Input.Applied_Force_Data.period = 2*pi/Nc_Inputs.frequency;
                            Eom_Input.Applied_Force_Data.amplitude = Nc_Inputs.amplitude;
                            Eom_Input.Applied_Force_Data.type = Nc_Inputs.force_type;

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
            dir_name = obj.data_path + "rom_data_" + obj.id;
            if ~isfolder(dir_name)
                mkdir(dir_name);
            end
            save(rom_data,"Eom_Input");
            pre_dynamic_time = toc(pre_dynamic_time_start);
            log_message = sprintf("EoM precomputations: %.1f seconds" ,pre_dynamic_time);
            logger(log_message,3)
        end
        %-----------------------------------------------------------------%
        function L_eigenvectors = get_current_L_eigenvectors(obj)
            L_modes = obj.Dynamic_Validation_Data.current_L_modes;
            all_L_modes = obj.Model.low_frequency_modes;
            [~,L_map] = ismember(L_modes,all_L_modes);
            L_evec_all = load_data(obj.Model.low_frequency_eigenvectors);
            L_eigenvectors = L_evec_all(:,L_map);
        end
        %-----------------------------------------------------------------%
        function [eom,eom_dz,eom_dt] = get_equation_of_motion(obj,varargin)
            %-------------------------------------------------------------------------%
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            Damping_Data = [];
            Force_Data = [];

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "damping"
                        Damping_Data = keyword_values{arg_counter};
                    case "forcing"
                        Force_Data = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------%
            conservative = (isempty(Damping_Data) && isempty(Force_Data));

            if conservative
                Eom_Input = obj.get_solver_inputs("coco_backbone");
                eom = @(t,z) coco_eom(0,z,0,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
                eom_dz = @(t,z) coco_eom_dx(0,z,0,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
            else
                if isempty(Force_Data)
                    Force_Data.type = "none";
                end


                switch Damping_Data.damping_type
                    case "rayleigh"
                        damping = get_rayleigh_damping_matrix(Damping_Data,obj.Model);
                        Nonconservative_Input.damping = damping;
                end

                switch Force_Data.type
                    case "modal"
                        mode_map = F_Data.mode_number == obj.Model.reduced_modes;
                        Nonconservative_Input.mode_map = mode_map;
                    case "point force"
                        num_dofs = obj.Model.num_dof;
                        dof_map = zeros(num_dofs,1);
                        dof_map(obj.Model.node_mapping(:,1) == F_Data.dof) = 1;
                        Nonconservative_Input.amplitude_shape = dof_map;
                    case "none"
                        Force_Data.amplitude = 0;
                        Force_Data.frequency = 0;
                end
                Nonconservative_Input.amplitude = Force_Data.amplitude;
                Nonconservative_Input.frequency = Force_Data.frequency;
                Nonconservative_Input.force_type = Force_Data.type;

   

                Eom_Input = obj.get_solver_inputs("coco_frf",Nonconservative_Input);


                T = Eom_Input.Applied_Force_Data.period;
                amp = Eom_Input.Applied_Force_Data.amplitude;
                eom = @(t,z) coco_forced_eom(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);
                eom_dz = @(t,z) coco_forced_eom_dx(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);
                eom_dt = @(t,z) coco_forced_eom_dt(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);


            end
        end
        %-----------------------------------------------------------------%
        function low_frequency_displacement = create_modal_validation_polynomial(obj)
            L_eigenvectors = get_current_L_eigenvectors(obj);
            r_eigenvectors = obj.Model.reduced_eigenvectors;
            h_eigenvectors = [r_eigenvectors,L_eigenvectors];

            modal_transform = h_eigenvectors'*obj.Model.mass;
            L_Modes_Poly = modal_transform*obj.Physical_Displacement_Polynomial;

            % L_Modes_Validation_Poly = modal_transform*obj.Low_Frequency_Coupling_Gradient_Polynomial;
            % low_frequency_displacement = @(r,h) L_Modes_Poly.evaluate_polynomial(r) + L_Modes_Validation_Poly.evaluate_polynomial(r)*h;
            low_frequency_displacement = @(r,h) L_Modes_Poly.evaluate_polynomial(r) + h;
        end
        %-----------------------------------------------------------------%
    end

end