classdef Polynomial
    %Tensors of multidimensional polynomials 
    properties
        polynomial_degree
        coefficients
        input_limit

        num_element_coefficients
        num_independent_element_coefficients
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
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            
            %---------------------
            if coupling_type == "stiffness"
                
                output_shape = size(output_data);
                num_rows = output_shape(1);
                for iRow = 1:num_rows
                    for iCol = (iRow+1):num_rows
                        element_1 = squeeze(output_data(iRow,iCol,:));
                        element_2 = squeeze(output_data(iCol,iRow,:));
                        mean_element = (element_1 + element_2)/2;
                        output_data(iRow,iCol,:) = mean_element;
                        output_data(iCol,iRow,:) = mean_element;
                    end
                end
            end
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

            obj.input_limit = Polynomial.find_limits(input_data);

            
            %---------------------
            if minimum_output ~= 0
                output_range = max(output_data,[],2) - min(output_data,[],2);
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
                end
                obj.modeled_outputs = ~zero_index;
            end
            %---------------------
            [scale_factor,shift_factor] = Polynomial.transform_data(input_data,scale,shift);
            obj.scaling_factor = scale_factor;
            obj.shifting_factor = shift_factor;
            
            %---------------------
            Constraint = obj.parse_constraint(constraint_type);
            
            %---------------------
            switch coupling_type
                case {"stiffness","none"}
                    [coeffs,num_unconstrained_terms] = obj.fit_polynomial(input_data,output_data,Constraint);
                case "force"
                    Coupling = obj.get_force_coupling;
                    [coeffs,num_unconstrained_terms] = obj.fit_coupled_polynomial(input_data,output_data,Constraint,Coupling);
            end
            %---------------------
            if output_reshaped
                obj.output_dimension = output_shape(1:2);
                coeffs = reshape(coeffs,[size(coeffs,1),output_shape(1:2)]);
            end


            %---------------------
            obj.coefficients = coeffs;
            obj.num_independent_element_coefficients = num_unconstrained_terms;
        end
        %-----------------------------------------------------------------%
        function [coeffs,num_unconstrained_terms] = fit_polynomial(obj,input_data,output_data,Constraint)
            %no coefficient coupling
            scale_factor = obj.scaling_factor;
            shift_factor = obj.shifting_factor;

            %--------------
            transformed_data = scale_factor.*(input_data - shift_factor); 
            
            %--------------
            [unconstrained_input_matrix,output_correction] = obj.apply_constraints(transformed_data,Constraint);

            %--------------
            output_data = output_data' + output_correction;

            %-------------
            coeffs = lsqminnorm(unconstrained_input_matrix,output_data,'warn');
            num_unconstrained_terms = size(coeffs,1);
            %--------------
            if ~isempty(Constraint.terms)
                coeffs = obj.reconstruct_coefficients(coeffs,Constraint);
            end
        end
        %-----------------------------------------------------------------%
        function [coeffs,num_unconstrained_terms] = fit_coupled_polynomial(obj,input_data,output_data,Constraint,Coupling)
            %conservative coefficient coupling
            input_index = obj.input_order;
            scale_factor = obj.scaling_factor;
            shift_factor = obj.shifting_factor;

            %--------------
            transformed_data = scale_factor.*(input_data - shift_factor); 
            
            %--------------
            [unconstrained_input_matrix,output_correction] = obj.apply_constraints(transformed_data,Constraint);
            
            %--------------
            output_data = output_data' + output_correction;


            %--------------
            %coupling
            coupling_pattern = Coupling.coupling_pattern;
            coupling_scale_factors = Coupling.coupling_scale_factors;
            num_inputs = size(input_index,2);
            if ~isempty(Constraint.terms)
                constrained_terms = Constraint.terms;
                coupling_pattern(constrained_terms,:) = [];
                coupling_scale_factors(constrained_terms,:) = [];
            end

            coupling_pattern = Polynomial.relabel_coupling_pattern(coupling_pattern);
            num_unique_coeffs = max(coupling_pattern,[],"all");
            num_points = size(unconstrained_input_matrix,1);
            coupled_input_matrix = zeros(num_inputs*num_points,num_unique_coeffs);

            combination_scaling = zeros(num_unique_coeffs,num_inputs);
            new_coupling_scale_factors = zeros(size(coupling_scale_factors));
            for iCoeff = 1:num_unique_coeffs
                [term,input] = find(coupling_pattern == iCoeff);
                input_col = zeros(num_points,num_inputs);
                col_scale = zeros(2,num_inputs);
                num_coupled_terms = length(input);
                for iInput = 1:num_coupled_terms

                    col_scale(1,input(iInput)) = coupling_scale_factors(term(iInput),input(iInput));
                    % col_scale(2,input(iInput)) = 1/coeff_scaling(term(iInput));
                    col_scale(2,input(iInput)) = scale_factor(input(iInput));

                    input_col(:,input(iInput)) = unconstrained_input_matrix(:,term(iInput));

                end
                
                if isscalar(input)
                    combination_scaling(iCoeff,input) = 1;
                    coupling_scale_factors(term(iInput),input(iInput)) = 1;
                    new_coupling_scale_factors(term(iInput),input(iInput)) = 1;
                else
                    col_prod = prod(col_scale);
                    combination_scaling(iCoeff,:) = col_prod./mean(col_prod);
                    for iInput = 1:num_coupled_terms
                        new_coupling_scale_factors(term(iInput),input(iInput)) = col_prod(input(iInput))./mean(col_prod);
                    end
                end
                scaled_input_col = input_col.*combination_scaling(iCoeff,:);
                coupled_input_matrix(:,iCoeff) = reshape(scaled_input_col,num_inputs*num_points,1);
            end
            
            coupled_output_data = reshape(output_data,numel(output_data),1);

           
            %-------------
            coupled_coeffs = lsqminnorm(coupled_input_matrix,coupled_output_data,'warn');
            num_unconstrained_terms = size(coupling_pattern,1);
            %--------------
            %Reconstruct coupling
            coeffs = coupled_coeffs(coupling_pattern).*new_coupling_scale_factors;
             
            % Reconstruct constraints
            if ~isempty(Constraint.terms)
                coeffs = obj.reconstruct_coefficients(coeffs,Constraint);
            end
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
            input_data = scale_factor.*(input_data-shift_factor);

            input_matrix = Polynomial.get_input_matrix(input_data,input_index);
            
            outputs = obj.output_dimension;
            if nargin == 2
                if isscalar(outputs)
                    output_data = input_matrix*coeffs;
                else
                    output_data = squeeze(tensorprod(input_matrix,coeffs,2,1));
                end

            else

                if isscalar(outputs)
                    output_data = input_matrix*coeffs(:,output);

                else
                    switch size(output,2)
                        case 1
                            coeff_size = size(coeffs);

                            coeffs = reshape(coeffs,coeff_size(1),prod(coeff_size(2:3)));
                            output_data = input_matrix*coeffs(:,output);
                        case 2
                            output_data = squeeze(tensorprod(input_matrix,coeffs(:,output(:,1),output(:,2)),2,1));
                    end
                end
            end

            
            if size(output_data,1) ~= obj.output_dimension(1)
                switch ndims(output_data)
                    case 2
                        output_data = output_data';
                    case 3
                      
                end
            end

        end
        %-----------------------------------------------------------------%
        function plot_polynomial(obj,ax,plotted_outputs)
            PLOT_RESOLUTION = 101;

            if ~exist("plotted_outputs","var")
                plotted_outputs = (1:obj.output_dimension)';
            end

            if exist("ax","var") && isempty(ax)
                clear ax
            end

            if exist("ax","var") && ~iscell(ax)
                ax_scalar = ax;
                ax = cell(1,1);
                ax{1,1} = ax_scalar;
            end


            
            if ismatrix(plotted_outputs)
                num_plotted_outputs = size(plotted_outputs,1);
            else
                num_plotted_outputs = length(plotted_outputs);
            end
            
            num_inputs = obj.input_dimension;
            limits = obj.input_limit;
            

            switch num_inputs
                case 1
                    x = linspace(limits(1),limits(2),PLOT_RESOLUTION);
                    

                    if ~exist("ax","var")
                        figure
                        tiledlayout(1,num_plotted_outputs)
                        for iOutput = 1:num_plotted_outputs
                            ax{1,iOutput} = nexttile;
                            xlabel("x")
                            if size(plotted_outputs,2) > 1 
                                ylabel("y_{(" + plotted_outputs(iOutput,1) + "," + plotted_outputs(iOutput,2) + ")}")
                            else
                                ylabel("y_" + plotted_outputs(iOutput))
                            end
                            box on
                        end
                    end

                    for iOutput = 1:num_plotted_outputs
                        if ismatrix(plotted_outputs)
                            plotted_output = plotted_outputs(iOutput,:);
                        else
                            plotted_output = plotted_outputs(iOutput);
                        end

                        y = obj.evaluate_polynomial(x,plotted_output);

                        hold(ax{1,iOutput},"on")
                        plot(ax{1,iOutput},x,y)
                        hold(ax{1,iOutput},"off")
                    end

                case 2
 
                    num_points = PLOT_RESOLUTION^2;

                    poly_bound = polyshape(limits');
                    [x_lim,y_lim] = boundingbox(poly_bound);


                    x = linspace(x_lim(1),x_lim(2),PLOT_RESOLUTION);
                    y = linspace(y_lim(1),y_lim(2),PLOT_RESOLUTION);
                    [X,Y] = meshgrid(x,y);

                    X_BC = nan(PLOT_RESOLUTION);
                    Y_BC = nan(PLOT_RESOLUTION);
                    for iCol = 1:PLOT_RESOLUTION
                        x_vec = [X(:,iCol),Y(:,iCol)];
                        valid_point = isinterior(poly_bound,x_vec);
                        X_BC(valid_point,iCol) = X(valid_point,iCol);
                        Y_BC(valid_point,iCol) = Y(valid_point,iCol);
                    end
                

                    if ~exist("ax","var")
                        figure
                        tiledlayout(1,num_plotted_outputs)
                        for iOutput = 1:num_plotted_outputs
                            ax{1,iOutput} = nexttile;
                            xlabel("x_1")
                            ylabel("x_2")
                            if size(plotted_outputs,2) > 1
                                zlabel("y_{(" + plotted_outputs(iOutput,1) + plotted_outputs(iOutput,2) + ")}")
                            else
                                zlabel("y_" + plotted_outputs(iOutput))
                            end
                            box on
                        end
                    end

                    for iOutput = 1:num_plotted_outputs
                        if ismatrix(plotted_outputs)
                            plotted_output = plotted_outputs(iOutput,:);
                        else
                            plotted_output = plotted_outputs(iOutput);
                        end

                        X_array = reshape(X_BC,1,num_points);
                        Y_array = reshape(Y_BC,1,num_points);
                        Z_array = obj.evaluate_polynomial([X_array;Y_array],plotted_output);
                        Z_BC = reshape(Z_array,size(X_BC));
                        
                        hold(ax{1,iOutput},"on")
                        mesh(ax{1,iOutput},X_BC,Y_BC,Z_BC)
                        hold(ax{1,iOutput},"off")
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

            int_coeffs = (obj.coefficients./squeeze(int_scale_factor)).*(1./obj.scaling_factor');

    
            
            coeffs_int = zeros(Poly_Int.num_element_coefficients,Poly_Int.output_dimension);
            for iTerm_int = 2:Poly_Int.num_element_coefficients
                term = Poly_Int.input_order(iTerm_int,:);
                % [test,term_index] = ismember(term,integrated_terms,'rows');
                matching_terms = cellfun(@(term_row) isequal(term_row,term),integrated_terms);
                coeff_index = reshape(matching_terms,size(int_coeffs));
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
            Poly_Int.coefficients(1,:) = -int_constant; 
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



            coeffs_size = num2cell([Poly_Diff.num_element_coefficients,Poly_Diff.output_dimension]);
            coeffs_diff = zeros(coeffs_size{:});
            coeffs_index = cell(size(coeffs_size));
            for iDim = 1:length(coeffs_size)
                coeffs_index{1,iDim} = 1:coeffs_size{1,iDim};
            end
            for iTerm_diff = 1:Poly_Diff.num_element_coefficients
                coeffs_index{1,1} = iTerm_diff;
                term = Poly_Diff.input_order(iTerm_diff,:);

                for iDiff = 1:Poly_Diff.input_dimension
                    coeffs_index{1,end} = iDiff;

                    [~,term_index] = ismember(term,diff_input_index(:,:,iDiff),'rows');

                    matched_coeffs = diff_coeffs(term_index,:,:)*diff_scale_factor(term_index,1,iDiff)*obj.scaling_factor(iDiff);
                    
                    coeffs_diff(coeffs_index{:}) = matched_coeffs;
                end

                % term_index = mod(term_index,int_input_size(1));
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
                    Constraint.values = constraint_type{1,2}';
                case "linear_disp"
                    Constraint.terms = 1:(num_inputs+1);
                    Constraint.values = zeros(num_inputs+1,num_outputs);
                    constraint_value = constraint_type{1,2};
                    
                    if isscalar(constraint_value)
                        constraint_value = ones(num_inputs,num_outputs)*constraint_value;
                    elseif size(constraint_value,1) == num_outputs
                        constraint_value = constraint_value';
                    end

                    for iMode = 1:num_inputs
                        Constraint.values(iMode+1,:) = constraint_value(iMode,:)/scale_factor(iMode);
                    end
                case "linear_force"
                    Constraint.terms = 1:(num_inputs+1);
                    constraint_values = constraint_type{1,2};

                    for iMode = 1:num_inputs
                        Constraint.values(1,iMode) = 0;
                        Constraint.values(1+iMode,iMode) = constraint_values(iMode)/scale_factor(iMode);
                    end

            end
        end
        %-----------------------------------------------------------------%
        function [unconstrained_input_matrix,output_correction] = apply_constraints(obj,transformed_data,Constraint)
            input_index = obj.input_order;
            num_inputs = obj.input_dimension;
            num_outputs = obj.output_dimension;
            scale_factor = obj.scaling_factor;
            shift_factor = obj.shifting_factor;
            
            num_terms = size(input_index,1);
            num_loadcases = size(transformed_data,2);

            terms = 1:num_terms;
            unconstrained_terms = terms;
            if ~isempty(Constraint.terms)
                constrained_terms = Constraint.terms;
                unconstrained_terms(constrained_terms) = [];
            else
                constrained_terms = [];
            end
            
            input_matrix = Polynomial.get_input_matrix(transformed_data,input_index);
            output_correction = zeros(num_loadcases,num_outputs);

            if isempty(Constraint.terms)
                unconstrained_input_matrix = input_matrix;
                return
            end

            
            constraint_values = Constraint.values;
            transformation_prod = -scale_factor.*shift_factor;
            % constant constraint
            constrained_input_matrix = input_matrix(:,1);
            unconstrained_input_matrix = input_matrix(:,2:end); 

            trans_prod_input_matrix = Polynomial.get_input_matrix(transformation_prod,input_index);
            unconstrained_scaling = -trans_prod_input_matrix(2:end);
            constrained_scaling = trans_prod_input_matrix(1);

            constrained_term_col = constrained_input_matrix(:,constrained_scaling == 1);
            unconstrained_input_matrix = unconstrained_input_matrix + constrained_term_col.*unconstrained_scaling;
            for iOutput = 1:num_outputs
                output_correction(:,iOutput) = output_correction(:,iOutput) - constrained_term_col*constraint_values(constrained_scaling==1,iOutput);
            end
            
            if isscalar(Constraint.terms)
                return
            end
            % linear constraint
            constrained_input_matrix = [constrained_input_matrix,unconstrained_input_matrix(:,1:num_inputs)];
            unconstrained_input_matrix = unconstrained_input_matrix(:,(num_inputs+1):end); 
            [diff_input_index,diff_scale_factor] = Polynomial.differentiate_input_index(input_index,1:num_inputs);
            for iDiff = 1:num_inputs
                di_input_index = diff_input_index(:,:,iDiff);
                diff_input_matrix = Polynomial.get_input_matrix(transformation_prod,di_input_index);
                diff_input_matrix(isnan(diff_input_matrix)) = 0;
                linear_constraint_scaling = diff_input_matrix.*(squeeze(diff_scale_factor(:,1,iDiff))');
                unconstrained_scaling = -linear_constraint_scaling(unconstrained_terms);
                constrained_scaling = linear_constraint_scaling(constrained_terms);

                constrained_term_col = constrained_input_matrix(:,constrained_scaling == 1);
                unconstrained_input_matrix = unconstrained_input_matrix + constrained_term_col.*unconstrained_scaling;
                for iOutput = 1:num_outputs
                    output_correction(:,iOutput) = output_correction(:,iOutput) - constrained_term_col*constraint_values(constrained_scaling==1,iOutput);
                end
            end
        end
        %-----------------------------------------------------------------%
        function coeffs = reconstruct_coefficients(obj,coeffs,Constraint)
            input_index = obj.input_order;
            num_terms = size(input_index,1);
            num_inputs = obj.input_dimension;
            num_outputs = obj.output_dimension;
            scale_factor = obj.scaling_factor;
            shift_factor = obj.shifting_factor;

            terms = 1:num_terms;
            unconstrained_terms = terms;

            constrained_terms = Constraint.terms;
            unconstrained_terms(constrained_terms) = [];


            transformation_prod = -scale_factor.*shift_factor;
            % linear constraint
            if ~isscalar(Constraint.terms)
                lin_coeffs = zeros(num_inputs,num_outputs);
                lin_constraint = Constraint.values(2:end,:);

                [diff_input_index,diff_scale_factor] = Polynomial.differentiate_input_index(input_index,1:num_inputs);
                for iDiff = 1:num_inputs
                    di_input_index = diff_input_index(:,:,iDiff);
                    diff_input_matrix = Polynomial.get_input_matrix(transformation_prod,di_input_index);
                    diff_input_matrix(isnan(diff_input_matrix)) = 0;
                    linear_constraint_scaling = diff_input_matrix.*(squeeze(diff_scale_factor(:,1,iDiff))');
                    unconstrained_scaling = -linear_constraint_scaling(unconstrained_terms);
                    lin_coeffs(iDiff,:) = lin_constraint(iDiff,:) + unconstrained_scaling*coeffs;
                end
                coeffs = [lin_coeffs;coeffs];
            end
            % constraint constraint
            
            const_constraint = Constraint.values(1,:);
            trans_prod_input_matrix = Polynomial.get_input_matrix(transformation_prod,input_index);
            unconstrained_scaling = -trans_prod_input_matrix(2:end);
            const_coeffs = const_constraint + unconstrained_scaling*coeffs;

            coeffs = [const_coeffs;coeffs];
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

            input_matrix = zeros(num_points,num_terms);
            for iTerm = 1:num_terms
                term_powers = input_index(iTerm,:);
                term_data = prod(input_data.^(term_powers),2);
                input_matrix(:,iTerm) = term_data;
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
        function limits = find_limits(input_data)
            num_inputs = size(input_data,1);
            switch num_inputs
                case 1
                    limits = [min(input_data),max(input_data)];
                case 2
        
                    bound_index = boundary(input_data(1,:)',input_data(2,:)',0.1);
                    limits = input_data(:,bound_index);
                otherwise
                    limits = [];
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
                shift_factor = mean(input_data,2);
            else
                shift_factor = zeros(num_inputs,1);
            end
            % transformed_data = scale_factor.*(input_data - shift_factor);
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
    end

end