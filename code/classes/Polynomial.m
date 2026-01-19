classdef Polynomial
    %Tensors of multidimensional polynomials 
    properties
        polynomial_degree
        coefficients
        input_limit

        num_element_coefficients
        num_independent_element_coefficients
        num_fitted_coefficients
        
        input_dimension
        output_dimension
        
        input_order

        scaling_factor
        shifting_factor

        modeled_outputs
    end
    methods
        function obj = Polynomial(input_data,output_data,degree,varargin)
             %Optional argumanents
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);
            
            coupling_type = "none";
            constraint_type = {"none",[]}; 
            scale = true;
            shift = true;
            minimum_output = 0;
            num_applied_forces = 0;

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "constraint"
                        constraint_type = keyword_values{arg_counter};
                    case "scale"
                        scale = keyword_values{arg_counter};
                    case "shift"
                        shift = keyword_values{arg_counter};
                    case "coupling"
                        coupling_type = keyword_values{arg_counter};
                    case "minimum_output"
                        minimum_output = keyword_values{arg_counter};
                    case "nonconservative"
                        num_applied_forces = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            
            %---------------------
            % if coupling_type == "stiffness"
            % 
            %     output_shape = size(output_data);
            %     num_rows = output_shape(1);
            %     for iRow = 1:num_rows
            %         for iCol = (iRow+1):num_rows
            %             element_1 = squeeze(output_data(iRow,iCol,:));
            %             element_2 = squeeze(output_data(iCol,iRow,:));
            %             mean_element = (element_1 + element_2)/2;
            %             output_data(iRow,iCol,:) = mean_element;
            %             output_data(iCol,iRow,:) = mean_element;
            %         end
            %     end
            % end
            output_reshaped = ndims(output_data) == 3;
            if output_reshaped
                output_shape = size(output_data);
                output_data = reshape(output_data,[prod(output_shape(1:2)),output_shape(3)]);
                if constraint_type{1,1} ~= "none"
                    constraint_type{1,2} = reshape(constraint_type{1,2},[prod(output_shape(1:2)),1]);
                end
            end
            %---------------------

            obj.polynomial_degree = degree;

            obj.input_dimension = size(input_data,1);
            obj.output_dimension = size(output_data,1);

            num_coeffs = Polynomial.input_combinations(degree,obj.input_dimension);
            input_index = Polynomial.get_input_index(degree,obj.input_dimension);

            obj.num_element_coefficients = num_coeffs;
            obj.input_order = input_index;



            
            %---------------------
            if minimum_output ~= 0
                output_range = max(abs(output_data),[],2);
                % max_output = max(abs(output_data),[],2);
                % zero_index = max_output < minimum_output;
                zero_index = output_range < minimum_output;
                
                switch constraint_type{1,1}
                    case "none"
                        output_data(zero_index,:) = repmat(mean(output_data(zero_index,:),2),1,size(output_data,2));
                    case "constant"
                        output_data(zero_index,:) = repmat(constraint_type{1,2}(zero_index,1),1,size(output_data,2));
                    case {"linear","linear_disp"}
                        output_data(zero_index,:) = 0;
                        constraint_value = constraint_type{1,2};
                        constraint_value(zero_index,:) = 0;
                        constraint_type{1,2} = constraint_value;
                end
                obj.modeled_outputs = ~zero_index;
            end
     
            %---------------------
            [scale_factor,shift_factor] = Polynomial.transform_data(input_data,scale,shift);
            obj.scaling_factor = scale_factor;
            obj.shifting_factor = shift_factor;
                       %---------------------
            Constraint = obj.parse_constraint(constraint_type);
            Constraint.coupling = coupling_type;
            obj.input_limit = Polynomial.find_limits(input_data,Constraint,num_applied_forces);
        
            
            %---------------------
            [coeffs,Num_Unconstrained_Terms] = obj.fit_polynomial(input_data,output_data,Constraint);
            switch coupling_type
                case "none"
                    obj.num_independent_element_coefficients = Num_Unconstrained_Terms.regression_size;
                    obj.num_fitted_coefficients = Num_Unconstrained_Terms.regression_size*size(coeffs,2);
                case "force"
                    obj.num_independent_element_coefficients = Num_Unconstrained_Terms.coupling_size;
                    obj.num_fitted_coefficients = Num_Unconstrained_Terms.regression_size;
                case "stiffness"
                    obj.num_independent_element_coefficients = Num_Unconstrained_Terms.regression_size;
                    obj.num_fitted_coefficients = Num_Unconstrained_Terms.regression_size*Num_Unconstrained_Terms.coupling_size;
            end
            %---------------------
            if output_reshaped
                obj.output_dimension = output_shape(1:2);
                coeffs = reshape(coeffs,[output_shape(1:2),size(coeffs,2)]);
            end


            %---------------------
            obj.coefficients = coeffs;
            
        end
        %-----------------------------------------------------------------%
        function [coeffs,Num_Unconstrained_Terms] = fit_polynomial(obj,input_data,output_data,Constraint)
            scale_factor = obj.scaling_factor;
            shift_factor = obj.shifting_factor;

            %--------------
            transformed_data = scale_factor.*(input_data + shift_factor); 
            
            %--------------
            [unconstrained_input_matrix,output_correction,Coeff_Data] = obj.apply_constraints(transformed_data,Constraint);
            output_data = output_data + output_correction;
            
            %-------------
            [output_data,unconstrained_input_matrix,Coupling_Data] = obj.apply_coupling(output_data,unconstrained_input_matrix,Constraint);

            %-------------
            coeffs = lsqminnorm(unconstrained_input_matrix',output_data','warn')';
            Num_Unconstrained_Terms.regression_size = size(coeffs,2);
            if ~isempty(Coupling_Data)
                Num_Unconstrained_Terms.coupling_size = Coupling_Data.couping_size;
            end
            %--------------
            coeffs = obj.reconstruct_coefficients(coeffs,Constraint,Coeff_Data,Coupling_Data);
        end
        %-----------------------------------------------------------------%
        
        %Using the polynomial
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function output_data = evaluate_polynomial(obj,input_data,output)
            coeffs = obj.coefficients;
            
            input_index = obj.input_order;

            
            shift_factor = obj.shifting_factor;
            scale_factor = obj.scaling_factor;
            input_data = scale_factor.*(input_data+shift_factor);

            input_matrix = Polynomial.get_input_matrix(input_data,input_index);
            
            outputs = obj.output_dimension;
            if nargin == 2
                if isscalar(outputs)
                    output_data = coeffs*input_matrix;
                else
                    output_data = squeeze(tensorprod(coeffs,input_matrix,ndims(coeffs),1));
                end

            else

                if isscalar(outputs)
                    output_data = coeffs(output,:)*input_matrix;

                else
                    switch size(output,2)
                        case 1
                            coeff_size = size(coeffs);

                            coeffs = reshape(coeffs,prod(coeff_size(1:2)),coeff_size(3));
                            output_data = coeffs(output,:)*input_matrix;
                        case 2
                            output_data = squeeze(tensorprod(coeffs(output(:,1),output(:,2),:),input_matrix,ndims(coeffs),1));
                    end
                end
            end

            
            % if size(output_data,1) ~= obj.output_dimension(1)
            %     switch ndims(output_data)
            %         case 2
            %             output_data = output_data';
            %         case 3
            % 
            %     end
            % end

        end
        %-----------------------------------------------------------------%
        function ax = plot_polynomial(obj,varargin)

            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            ax = [];
            tag = "";
            total_outputs = prod(obj.output_dimension);
            plotted_outputs = (1:total_outputs)';
            Potential_Poly = [];
            energy_limit = [];
            colour_num = 1;

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "axes"
                        ax = keyword_values{arg_counter};
                        if ~iscell(ax)
                            ax = {ax};
                        end
                    case "tag"
                        tag = keyword_values{arg_counter};
                    case "outputs"
                        plotted_outputs = keyword_values{arg_counter};
                    case "potential"
                        potential_data = keyword_values{arg_counter};
                        Potential_Poly = potential_data{1};
                        energy_limit = potential_data{2};
                    case {"color","colour"}
                        colour_num = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %------------


            PLOT_RESOLUTION = 101;




            num_plotted_outputs = size(plotted_outputs,1);


            num_inputs = obj.input_dimension;
            limits = obj.input_limit;


            switch num_inputs
                case 1
                    x = linspace(limits(1),limits(2),PLOT_RESOLUTION);


                    if isempty(ax)
                        figure
                        tiledlayout(1,num_plotted_outputs)
                        ax = cell(num_plotted_outputs,1);
                        for iOutput = 1:num_plotted_outputs
                            ax{iOutput} = nexttile;
                            xlabel("x")
                            if size(plotted_outputs,2) > 1 
                                ylabel("y_{(" + plotted_outputs(iOutput,1) + "," + plotted_outputs(iOutput,2) + ")}")
                            else
                                ylabel("y_{" + plotted_outputs(iOutput)+ "}")
                            end
                            box on
                        end
                    end

                    for iOutput = 1:num_plotted_outputs
                        hold_state = ishold(ax{iOutput});

                        plotted_output = plotted_outputs(iOutput,:);


                        y = obj.evaluate_polynomial(x,plotted_output);
                        
                        hold(ax{iOutput},"on")
                        plot(ax{iOutput},x,y,"Color",get_plot_colours(colour_num),"Tag",tag)
                        hold(ax{iOutput},"off")

                        if hold_state
                            hold(ax{iOutput},"on");
                        end
                    end

                case 2
                    [x_grid,y_grid] = get_2d_plotting_data(obj,Potential_Poly,energy_limit);
                    num_points = numel(x_grid);
                    
                    if isempty(ax)
                        figure
                        tiledlayout("flow")
                        ax = cell(num_plotted_outputs);
                        for iOutput = 1:num_plotted_outputs
                            ax{iOutput} = nexttile;
                            xlabel("x_1")
                            ylabel("x_2")
                            if size(plotted_outputs,2) > 1
                                zlabel("y_{(" + plotted_outputs(iOutput,1) + plotted_outputs(iOutput,2) + ")}")
                            else
                                zlabel("y_{" + plotted_outputs(iOutput) + "}")
                            end
                            box on
                        end
                    end

                    for iOutput = 1:num_plotted_outputs
                        hold_state = ishold(ax{iOutput});
                        plotted_output = plotted_outputs(iOutput,:);

                        x_lin = reshape(x_grid,1,num_points);
                        y_lin = reshape(y_grid,1,num_points);
                        z_lin = obj.evaluate_polynomial([x_lin;y_lin],plotted_output);
                        z_grid = reshape(z_lin',size(x_grid));
                        
                        hold(ax{iOutput},"on")
                        mesh(ax{iOutput},x_grid,y_grid,z_grid,"Tag",tag)
                        hold(ax{iOutput},"off")

                        if hold_state
                            hold(ax{iOutput},"on");
                        end
                    end
            end
            
        
        
        end
        %-----------------------------------------------------------------%
        
        
        %Polynomial Manipulation
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function Poly_Int = integrate_polynomial(obj)
            Poly_Int = obj;
            Poly_Int.polynomial_degree = obj.polynomial_degree + 1;
            Poly_Int.output_dimension = 1;
            
            input_index_int = Polynomial.get_input_index(obj.polynomial_degree+1,obj.input_dimension);
            Poly_Int.input_order = input_index_int;
            Poly_Int.num_element_coefficients = size(input_index_int,1);

            [int_input_index,int_scale_factor] = Polynomial.integrate_input_index(obj.input_order);
            % int_input_size = size(int_input_index);
            integrated_terms = reshape(permute(int_input_index,[1,3,2]),[size(int_input_index,1)*size(int_input_index,3),size(int_input_index,2)]);
            integrated_terms = num2cell(integrated_terms,2);
            
            int_scale_factor = squeeze(int_scale_factor)';
            int_coeffs = (obj.coefficients./int_scale_factor)./obj.scaling_factor;

    
            
            coeffs_int = zeros(Poly_Int.output_dimension,Poly_Int.num_element_coefficients);
            for iTerm_int = 2:Poly_Int.num_element_coefficients
                term = Poly_Int.input_order(iTerm_int,:);
                % [test,term_index] = ismember(term,integrated_terms,'rows');
                matching_terms = cellfun(@(term_row) isequal(term_row,term),integrated_terms);
                coeff_index = reshape(matching_terms,flip(size(int_coeffs)))';
                matched_coeffs = int_coeffs(coeff_index);

                if isscalar(uniquetol(matched_coeffs))
                    coeffs_int(iTerm_int) = matched_coeffs(1);
                else
                    error("Potential not consistent")
                end
                % term_index = mod(term_index,int_input_size(1));
            end

            Poly_Int.coefficients = coeffs_int;
            int_constant = Poly_Int.evaluate_polynomial(zeros(Poly_Int.input_dimension,1));
            Poly_Int.coefficients(:,1) = -int_constant; 
        end
        %-----------------------------------------------------------------%
        function Poly_Diff = differentiate_polynomial(obj)
            Poly_Diff = obj;
            Poly_Diff.polynomial_degree = obj.polynomial_degree - 1;
            Poly_Diff.output_dimension = [obj.output_dimension,obj.input_dimension];

            input_index_diff = Polynomial.get_input_index(Poly_Diff.polynomial_degree,Poly_Diff.input_dimension);
            Poly_Diff.input_order = input_index_diff;
            Poly_Diff.num_element_coefficients = size(input_index_diff,1);

            [diff_input_index,diff_scale_factor] = Polynomial.differentiate_input_index(obj.input_order,1:obj.input_dimension);
            % diff_input_size = size(diff_input_index);
            % differentiated_terms = reshape(permute(diff_input_index,[1,3,2]),[diff_input_size(1)*diff_input_size(3),diff_input_size(2)]);
            % differentiated_terms = num2cell(differentiated_terms,2);
            
            diff_coeffs = obj.coefficients;
            
            diff_index = size(Poly_Diff.output_dimension,2);


            coeffs_size = num2cell([Poly_Diff.output_dimension,Poly_Diff.num_element_coefficients]);
            coeffs_diff = zeros(coeffs_size{:});
            coeffs_index = cell(size(coeffs_size));
            for iDim = 1:length(coeffs_size)
                coeffs_index{1,iDim} = 1:coeffs_size{1,iDim};
            end
            for iTerm_diff = 1:Poly_Diff.num_element_coefficients
                coeffs_index{1,end} = iTerm_diff;
                term = Poly_Diff.input_order(iTerm_diff,:);

                for iDiff = 1:Poly_Diff.input_dimension
                    coeffs_index{1,diff_index} = iDiff;

                    [~,term_index] = ismember(term,diff_input_index(:,:,iDiff),'rows');
                    

                    matched_coeffs = diff_coeffs(:,term_index)*diff_scale_factor(term_index,1,iDiff)*obj.scaling_factor(iDiff);
                    
                    coeffs_diff(coeffs_index{:}) = matched_coeffs;
                end


            end
            
            Poly_Diff.coefficients = coeffs_diff;
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        %Overloading
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function sz = size(obj,dim)
            sz = obj.output_dimension;
            if isscalar(sz)
                sz = [sz,1];
            end
            if nargin == 2
                sz = sz(dim);
            end
        end
        %-----------------------------------------------------------------%
        function poly_two = mtimes(const,poly_one)
            poly_two = poly_one;
            coeffs = poly_two.coefficients;
            % new_coeffs = (const*coeffs')';
            new_coeffs = tensorprod(const,coeffs,ndims(const),1);
            poly_two.coefficients = new_coeffs;
            out_dim = size(new_coeffs);
            out_dim(end) = [];
            poly_two.output_dimension = out_dim;

            if ~isempty(poly_two.modeled_outputs)
                if size(out_dim,2) == 1
                    out_dim(2) = 1;
                end
                poly_two.modeled_outputs = true(out_dim);
            end
        end
        %-----------------------------------------------------------------%
        function obj = subpoly(obj,index)
            %overloading was problematic with function calls
            if isstring(index) && index == "all"

            else
                obj.coefficients = obj.coefficients(index,:,:);
                obj.output_dimension(1) = length(index);
                if ~isempty(obj.modeled_outputs)
                    obj.modeled_outputs = obj.modeled_outputs(index);
                end
                obj.num_fitted_coefficients = obj.num_independent_element_coefficients * obj.output_dimension;
            end
        end
        % function poly_sub = subsref(obj,index)
        % 
        % end
        % function poly_sub = subsindex(obj,index)
        % 
        % end
        %-----------------------------------------------------------------%

        %Helpers
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function Constraint = parse_constraint(obj,constraint_type)
     
            input_index = obj.input_order;
            num_inputs = size(input_index,2);
            num_outputs = obj.output_dimension;

            scale_factor = obj.scaling_factor;

            switch constraint_type{1,1}
                case "none"
                    Constraint.terms = [];
                case "constant"
                    Constraint.terms = 1;
                    Constraint.values = constraint_type{1,2};
                case "linear"
                    Constraint.terms = 1:(num_inputs+1);
                    Constraint.values = constraint_type{1,2};
                    if isscalar(Constraint.values)
                        Constraint.values = ones(obj.output_dimension,num_inputs+1)*Constraint.values;
                    end
                case "linear_disp"
                    Constraint.terms = 1:(num_inputs+1);
                    Constraint.values = zeros(num_outputs,num_inputs+1);
                    constraint_value = constraint_type{1,2};
                    
                    if isscalar(constraint_value)
                        constraint_value = ones(num_outputs,num_inputs)*constraint_value;
                    end

                    for iMode = 1:num_inputs
                        Constraint.values(:,iMode+1) = constraint_value(:,iMode)/scale_factor(iMode);
                    end
                case "linear_force"
                    Constraint.terms = 1:(num_inputs+1);
                    constraint_values = constraint_type{1,2};

                    for iMode = 1:num_inputs
                        Constraint.values(iMode,1) = 0;
                        Constraint.values(iMode,1+iMode) = constraint_values(iMode)/scale_factor(iMode);
                    end

            end
        end
        %-----------------------------------------------------------------%
        function [unconstrained_input_matrix,output_correction,Coeff_Data] = apply_constraints(obj,transformed_data,Constraint)
            input_index = obj.input_order;
            scale_factor = obj.scaling_factor;
            shift_factor = obj.shifting_factor;
            num_loadcases = size(transformed_data,2);
            num_outputs = obj.output_dimension;
            num_inputs = obj.input_dimension;
            num_terms = size(input_index,1);
     
            input_matrix = Polynomial.get_input_matrix(transformed_data,input_index);
            
            if isempty(Constraint.terms)
                output_correction = zeros(num_outputs,num_loadcases);
                unconstrained_input_matrix = input_matrix;
                Coeff_Data = [];
                return
            end
            % a_constrained = A + a_unconstrained * B
            constrained_terms = Constraint.terms;
            num_constrained_terms = size(constrained_terms,2);
            
            unconstrained_terms = 1:num_terms;
            unconstrained_terms(constrained_terms) = [];


            transformation_prod = scale_factor.*shift_factor;
            trans_prod_input = Polynomial.get_input_matrix(transformation_prod,input_index);
            switch num_constrained_terms
                case 1
                    A = Constraint.values;
                    B = - trans_prod_input(unconstrained_terms);
                    
                    output_correction = - A;
                    unconstrained_input_matrix = input_matrix(unconstrained_terms,:) + B;

                case (1+num_inputs)
                    trans_prod_diff_input = Polynomial.get_diff_input_matrix(transformation_prod,input_index);
                    trans_prod_input = [trans_prod_input,trans_prod_diff_input];
                    
                    constrained_trans_prod_input = trans_prod_input(constrained_terms,:);
                    A = Constraint.values/constrained_trans_prod_input;

                    unconstrained_trans_prod_input = trans_prod_input(unconstrained_terms,:);
                    B = - unconstrained_trans_prod_input/constrained_trans_prod_input;

                    output_correction = -A *input_matrix(constrained_terms,:);
                    unconstrained_input_matrix = input_matrix(unconstrained_terms,:) + B*input_matrix(constrained_terms,:);
            end
            Coeff_Data.constrained_coefficient = A;
            Coeff_Data.unconstrained_coefficient = B;

        end
        %-----------------------------------------------------------------%
        function [coupled_output_data,coupled_input_matrix,Coupling_Data] = apply_coupling(obj,output_data,unconstrained_input_matrix,Constraint)
            switch Constraint.coupling
                case "none"
                    coupled_output_data = output_data;
                    coupled_input_matrix = unconstrained_input_matrix;
                    Coupling_Data = [];
                case "force"
                    scale_factor = obj.scaling_factor;

                    constrained_terms = Constraint.terms;
                    num_constrained_terms = size(constrained_terms,2);
                    num_unconstrained_terms = size(unconstrained_input_matrix,1);
                    num_terms = num_unconstrained_terms + num_constrained_terms;

                    unconstrained_terms = 1:num_terms;
                    unconstrained_terms(constrained_terms) = [];

                    num_outputs = obj.output_dimension;
                    num_inputs = obj.input_dimension;

                    num_potential_coeffs = Polynomial.input_combinations(obj.polynomial_degree+1,num_inputs);
                    switch num_constrained_terms
                        case 1
                            num_potential_coeffs = num_potential_coeffs - Polynomial.input_combinations(1,num_inputs);
                        case (1 + num_inputs)
                            num_potential_coeffs = num_potential_coeffs - Polynomial.input_combinations(2,num_inputs);
                    end
                    %----

                    Coupling = obj.get_force_coupling;
                    

                    coupling_pattern = Coupling.coupling_pattern(unconstrained_terms,:);
                    coupling_pattern = coupling_pattern - min(coupling_pattern,[],"all") + 1;
                    
                    coupling_scale_factors = Coupling.coupling_scale_factors(unconstrained_terms,:);
                    
                    num_loadcases = size(output_data,2);
                    coupled_input_matrix = zeros(num_potential_coeffs,num_loadcases*num_outputs);
                    
                    coupling_col = zeros(1,num_outputs*num_unconstrained_terms);
                    coupling_row = zeros(1,num_outputs*num_unconstrained_terms);
                    coupling_value = zeros(1,num_outputs*num_unconstrained_terms);
                    for iOutput = 1:num_outputs
                        coupling_col_i = 1:num_unconstrained_terms;
                        coupling_row_i = coupling_pattern(coupling_col_i,iOutput);
                        coupling_value_i = coupling_scale_factors(coupling_col_i ,iOutput)*scale_factor(iOutput);
                        C_i = sparse(coupling_row_i,coupling_col_i,coupling_value_i,num_potential_coeffs,num_unconstrained_terms);

                        coupling_span = (num_loadcases*(iOutput-1) + 1):(num_loadcases*iOutput);
                        coupled_input_matrix(:,coupling_span) = C_i*unconstrained_input_matrix;

                        sparse_span = (num_unconstrained_terms*(iOutput-1) +1):(num_unconstrained_terms*iOutput);
                        coupling_col(sparse_span) = sparse_span;
                        coupling_row(sparse_span) = coupling_row_i;
                        coupling_value(sparse_span) = coupling_value_i;
                    end
                    Coupling_Data.C = sparse(coupling_row,coupling_col,coupling_value,num_potential_coeffs,num_unconstrained_terms*num_outputs);
                    Coupling_Data.couping_size = size(coupling_pattern,1);
                    %----
                    coupled_output_data = reshape(output_data',1,num_loadcases*num_outputs);
                case "stiffness"
                    num_terms = size(output_data,1);
                    matrix_dim = sqrt(num_terms);
                    num_unique_terms = matrix_dim*(matrix_dim+1)/2;

                    coupling = zeros(matrix_dim);
                    term_counter = 0;
                    for iRow = 1:matrix_dim
                        for iCol = iRow:matrix_dim
                            term_counter = term_counter + 1;
                            coupling(iRow,iCol) = term_counter;
                            coupling(iCol,iRow) = term_counter;
                        end
                    end
                    coupling_vec = reshape(coupling',num_terms,1);
                    
                    coupling_row = 1:num_terms;
                    coupling_col = coupling_vec(coupling_row);
                    coupling_value = ones(1,num_terms);

                    Coupling_Data.C = sparse(coupling_row,coupling_col,coupling_value,num_terms,num_unique_terms);
                    Coupling_Data.couping_size = num_unique_terms;
                    %---
                    coupled_input_matrix = unconstrained_input_matrix;
                    
                    %--
                    %check output data symmetry
                    num_loadcases = size(output_data,2);
                    coupled_output_data = zeros(num_unique_terms,num_loadcases);
                    for iTerm = 1:num_unique_terms
                        term_index = coupling_vec == iTerm;
                        coupled_data = output_data(term_index,:);
                        if nnz(term_index) > 1
                            is_symmetric = isapprox(coupled_data(1,:),coupled_data(2,:),"veryloose");
                            if any(~is_symmetric)
                                % coupled_data(:,~is_symmetric) = 0;
                                % asymmetric_terms = coupled_data(:,~is_symmetric);
                                if nnz(~is_symmetric) > num_loadcases/2
                                    coupled_data(:,:) = 0;
                                end
                            end
                            coupled_data = mean(coupled_data,1);
                        end

                        coupled_output_data(iTerm,:) = coupled_data;
                    end
            end

        end
        %-----------------------------------------------------------------%
        function coeffs = reconstruct_coefficients(obj,fitted_coeffs,Constraint,Coeff_Data,Coupling_Data)
            if isempty(Constraint.terms)
                coeffs = fitted_coeffs;
                return
            end

            switch Constraint.coupling
                case "force"
                    coupled_coeffs = fitted_coeffs;
                    fitted_coeffs = coupled_coeffs*Coupling_Data.C;
                    fitted_coeffs = reshape(fitted_coeffs,[],obj.output_dimension)';
                case "stiffness"
                    coupled_coeffs = fitted_coeffs;
                    fitted_coeffs = Coupling_Data.C*coupled_coeffs;
            end

            constraind_coeffs = Coeff_Data.constrained_coefficient + fitted_coeffs*Coeff_Data.unconstrained_coefficient;
            
            num_constrained_coeffs = size(constraind_coeffs,2);
            num_fitted_coeffs = size(fitted_coeffs,2);
            num_terms = num_fitted_coeffs + num_constrained_coeffs;
            num_outputs = size(fitted_coeffs,1);
            coeffs = zeros(num_outputs,num_terms);

            terms = 1:num_terms;
            unconstrained_terms = terms;
            unconstrained_terms(Constraint.terms) = [];
            coeffs(:,Constraint.terms) = constraind_coeffs;
            coeffs(:,unconstrained_terms) = fitted_coeffs;

        end
        %-----------------------------------------------------------------%
        function Coupling = get_force_coupling(obj)
            input_index = obj.input_order;
            num_inputs = size(input_index,2);

            [int_input_pattern,int_scale_factor] = Polynomial.integrate_input_index(input_index);

            int_input_index = Polynomial.get_input_index(obj.polynomial_degree+1,num_inputs);
            num_int_coefficients = size(int_input_index,1);

            num_terms = size(input_index,1);
            coupling_pattern = zeros(num_terms,num_inputs);
            coupling_scale_factors = zeros(num_terms,num_inputs);
            for iTerm = 1:num_int_coefficients
                term = int_input_index(iTerm,:);
                %check where this term exist in the derivatives
                for iInt = 1:num_inputs
                    term_index = ismember(int_input_pattern(:,:,iInt),term,'rows');

                    % coupling_pattern(iTerm,iDiff) = term_scale_factor;
                    coupling_pattern(term_index,iInt) = iTerm;
                    coupling_scale_factors(term_index,iInt) = int_scale_factor(term_index,1,iInt);
                end
            end

            Coupling.coupling_pattern = coupling_pattern;
            Coupling.coupling_scale_factors = coupling_scale_factors;
        end
        %-----------------------------------------------------------------%
        function Diff_Data = get_diff_data(obj,num_diff)
            starting_input_index = obj.input_order;
            data_scale_factor = obj.scaling_factor;


            scale_factor_size = size(starting_input_index);
            scale_factor_size(2) = 1;
            scale_factor = ones(scale_factor_size);
    
            num_inputs = obj.input_dimension;
            diff_inputs = 1:num_inputs;
            


            diff_input_index = cell(1,num_diff);
            diff_scale_factor = cell(1,num_diff);
            diff_mapping = cell(1,num_diff);
            
            input_index = starting_input_index;
            for iDiff = 1:num_diff
                [di_input_index,di_scale_factor] = Polynomial.differentiate_input_index(input_index,diff_inputs);
                
                
               
                scale_factor_dims = ndims(scale_factor);
                permutation_order = [1:(scale_factor_dims-1),scale_factor_dims+1,scale_factor_dims];
                pre_scale_factor = permute(scale_factor,permutation_order);
                scale_factor = pre_scale_factor.*di_scale_factor;

                scale_factor_index = cell(1,scale_factor_dims+1);
                % scale_factor_size = size(scale_factor);
                for iDim = 1:(scale_factor_dims+1)
                    scale_factor_index{1,iDim} = 1:size(scale_factor,iDim);
                end
                for iDiff_input = diff_inputs
                    scale_factor_index{1,max(3,scale_factor_dims)} = iDiff_input; %may not work for third derivative
                    scale_factor(scale_factor_index{:}) = scale_factor(scale_factor_index{:}) * data_scale_factor(iDiff_input);
                end
 
        

                % scale_factor_size = size(scale_factor);
                % pre_scale_factor = reshape(scale_factor,[scale_factor_size(1:(end-1)),1,scale_factor_size(end)]);
                % scale_factor = pre_scale_factor.*di_scale_factor;
                
                input_index = di_input_index;
                scale_factor(isnan(scale_factor)) = 0;

                diff_input_index{1,iDiff} = input_index;
                diff_scale_factor{1,iDiff} = squeeze(scale_factor);
                diff_mapping{1,iDiff} = Polynomial.get_input_index_mapping(starting_input_index,input_index);
            end
            Diff_Data.diff_input_index = diff_input_index;
            Diff_Data.diff_scale_factor = diff_scale_factor;
            Diff_Data.diff_mapping = diff_mapping;
        end
        %-----------------------------------------------------------------%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static)
        function num_coefficients = input_combinations(degree,num_inputs)
            num_coefficients = 1;
            for iTerm = 1:degree
                num_coefficients = num_coefficients + nchoosek(iTerm+num_inputs-1,num_inputs-1);
            end
        end
        %-----------------------------------------------------------------%
        function input_index = get_input_index(degree,num_inputs)
            num_coefficients = Polynomial.input_combinations(degree,num_inputs);
        
            %Power raised of each input for each term
            input_index = zeros(num_coefficients,num_inputs);
            term_counter = 1;
            for iDegree = 1:degree
                term_input_indices = nchoosek(1:num_inputs+iDegree-1,iDegree) - (0:iDegree-1);
                num_terms = size(term_input_indices,1);

                for iTerm = 1:num_terms
                    term_index = term_input_indices(iTerm,:);
                    term_counter = term_counter + 1;

                    for iInput = 1:num_inputs
                        input_index(term_counter,iInput) = sum(term_index == iInput);
                    end
                end
            end
        end
        %-----------------------------------------------------------------%
        function input_matrix = get_input_matrix(input_data,input_index)
            if size(input_index,2) ~= size(input_data,2)
                input_data = input_data';
            end

            num_points = size(input_data,1);
            num_terms = size(input_index,1);

            input_matrix = zeros(num_terms,num_points);
            for iTerm = 1:num_terms
                term_powers = input_index(iTerm,:);
                term_data = prod(input_data.^(term_powers),2);
                input_matrix(iTerm,:) = term_data';
            end
        end
        %-----------------------------------------------------------------%
        function diff_input_matrix = get_diff_input_matrix(input_data,input_index)
            num_inputs = size(input_index,2);
            num_terms = size(input_index,1);
            [diff_input_index,diff_scale_factor] = Polynomial.differentiate_input_index(input_index,1:num_inputs);
            diff_input_matrix = zeros(num_terms,num_inputs);
            for iDiff = 1:num_inputs
                di_input_index = diff_input_index(:,:,iDiff);
                unscaled_diff_input_matrix = Polynomial.get_input_matrix(input_data,di_input_index);
                unscaled_diff_input_matrix(isnan(unscaled_diff_input_matrix)) = 0;
                diff_input_matrix(:,iDiff) = unscaled_diff_input_matrix.*(squeeze(diff_scale_factor(:,1,iDiff)));
            end
        end
        %-----------------------------------------------------------------%
        function input_index_mapping = get_input_index_mapping(input_index_one,input_index_two)
            num_terms = size(input_index_one,1);
            % input_index_two_size = size(input_index_two);
            % mapping_output = input_index_two_size(3:end);
            mapping_output = size(input_index_two,3:max(ndims(input_index_two),3));
            num_rows = num_terms*prod(mapping_output);
            row_mapping = zeros(num_rows,1);
            
            permutating_index = [1,3:ndims(input_index_two),2];
            input_index_two_temp = permute(input_index_two,permutating_index);
            input_index_two_terms = reshape(input_index_two_temp,num_rows,[]);
            
            for iTerm = 1:num_terms
                term = input_index_one(iTerm,:);
                term_index = ismember(input_index_two_terms,term,'rows');
                row_mapping(term_index,1) = iTerm;
            end
            row_mapping(row_mapping == 0) = 1; %for the sake of efficient indexing (scale factor will take it to zero)
            input_index_mapping = reshape(row_mapping,[num_terms,mapping_output]);
            
        end
        %-----------------------------------------------------------------%
        function limits = find_limits(input_data,Constraint,num_applied_forces)

            if ismember(1,Constraint.terms)
                input_data = [input_data,zeros(size(input_data,1),1)];
            end
            num_inputs = size(input_data,1);
            switch num_inputs
                case 1
                    limits = [min(input_data),max(input_data)];
                case 2
        
                    bound_index = boundary(input_data(1,:)',input_data(2,:)',0.1);
                    limits = input_data(:,bound_index);
                    if num_applied_forces == 1 && false
                        applied_force_limits = limits(2,:);
                        limit_sign = sign(applied_force_limits);
                        applied_force_limits(limit_sign == 1) = pi/2;
                        applied_force_limits(limit_sign == -1) = -pi/2;
                        limits(2,:) = applied_force_limits;
                    end
                otherwise
                    limits = [min(input_data,[],2),max(input_data,[],2)];
            end
        end
        %-----------------------------------------------------------------%
        function [scale_factor,shift_factor] = transform_data(input_data,scale,shift)
            num_inputs = size(input_data,1);
            if scale
                scale_factor = 1./std(input_data,1,2);
            else
                scale_factor = ones(num_inputs,1);
            end
            if shift
                shift_factor = -mean(input_data,2);
            else
                shift_factor = zeros(num_inputs,1);
            end
            % transformed_data = scale_factor.*(input_data + shift_factor);
        end
        %-----------------------------------------------------------------%
        function [diff_input_index,diff_scale_factor] = differentiate_input_index(input_index,diff_inputs)
            num_diff_input = length(diff_inputs);
            num_inputs = size(input_index,2);

            input_index_size = [size(input_index),num_diff_input];
            diff_input_index = zeros(input_index_size);
            
            scale_factor_size = input_index_size;
            scale_factor_size(2) = 1;
            diff_scale_factor = zeros(scale_factor_size);

            num_dims = length(input_index_size);
            % num_dims = ndims(input_index_size);

            for iDiff = 1:num_diff_input
                diff_input = diff_inputs(iDiff);
                diff_index = diff_input == (1:num_inputs);
                
                scale_factor_span = cell(1,num_dims);
                diff_input_index_span = cell(1,num_dims);
                input_index_span = cell(1,num_dims-1);
                for iDim = 1:(num_dims-1)
                    scale_factor_span{1,iDim} = 1:scale_factor_size(iDim);
                    diff_input_index_span{1,iDim} = 1:input_index_size(iDim);
                    input_index_span{1,iDim} = 1:input_index_size(iDim);
                end
                input_index_span{1,end} = diff_index;
                scale_factor_span{1,end} = iDiff;
                diff_input_index_span{1,end} = iDiff;
                
                diff_scale_factor(scale_factor_span{:}) = input_index(input_index_span{:});
                diff_input_index(diff_input_index_span{:}) = input_index;

                diff_input_index_span{1,2} = diff_index;
                diff_input_index(diff_input_index_span{:}) = diff_input_index(diff_input_index_span{:}) - 1;

                zero_terms = diff_input_index(diff_input_index_span{:}) < 0;
                zero_index = false(input_index_size);
                
                for iInput = 1:num_inputs
                    diff_input_index_span{1,2} = iInput;
                    zero_index(diff_input_index_span{:}) = zero_terms;
                end
                % diff_input_index_span{1,1} = zero_terms(:,1,1);
               
                diff_input_index(zero_index) = NaN;

                % diff_scale_factor(:,1,iDiff) = input_index(:,diff_index);
                % diff_input_index(:,:,iDiff) = input_index;
                % diff_input_index(:,diff_index,iDiff) = diff_input_index(:,diff_index,iDiff) - 1;
                % 
                % zero_terms = diff_input_index(:,diff_index,iDiff) < 0;
                % diff_input_index(zero_terms,:,iDiff) = NaN;
            end
        end
        %-----------------------------------------------------------------%
        function [int_input_index,int_scale_factor] = integrate_input_index(input_index)
            num_inputs = size(input_index,2);
            int_input_index = zeros([size(input_index),num_inputs]);
            int_scale_factor = zeros([size(input_index,1),1,num_inputs]);

            for iInt = 1:num_inputs

                int_scale_factor(:,1,iInt) = input_index(:,iInt)+1;
                int_input_index(:,:,iInt) = input_index;
                int_input_index(:,iInt,iInt) = int_input_index(:,iInt,iInt) + 1;
            end
        end
        %-----------------------------------------------------------------%
        function coupling_pattern = relabel_coupling_pattern(coupling_pattern)
            
            order = reshape(coupling_pattern,numel(coupling_pattern),1);
            new_order = zeros(size(order));
            coefficient_counter = 0;
            for iOrder = 1:length(order)
                terms_index = order(iOrder) == order;
                if sum(new_order(terms_index)) == 0
                    coefficient_counter = coefficient_counter + 1;
                    new_order(terms_index) = coefficient_counter;
                   
                end
            end
            coupling_pattern = reshape(new_order,size(coupling_pattern));
        end
        %-----------------------------------------------------------------%
        
        
    end

end