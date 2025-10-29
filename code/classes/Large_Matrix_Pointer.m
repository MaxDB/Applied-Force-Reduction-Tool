classdef Large_Matrix_Pointer
    properties
        variable_path
        variable_name
        type
        matrix_size
    end

    methods
        function obj = Large_Matrix_Pointer(matrix,path,variable_name,varargin)
            %-------------------------------------------------------------------------
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);

            save = 1;

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "save"
                        save = keyword_values{arg_counter};
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %-------------------------------------------------------------------------
            if issparse(matrix)
                obj.type = "sparse";
            else
                obj.type = "full";
            end
            obj.matrix_size = size(matrix);

            obj.variable_path = path;
            obj.variable_name = variable_name;

            if save
                obj.save(matrix);
            end

            
        end
        %----
        function save(obj,matrix)
            file_path = obj.variable_path;
            if ~isfolder(file_path)
                mkdir(file_path)
            end
            file_path = file_path + obj.variable_name;
            save(file_path,"matrix")

        end
        %----
        function matrix = load(obj)
            if isfolder(obj.variable_path)
                load(obj.variable_path + obj.variable_name,"matrix");
            elseif isfile(obj.variable_path + ".mat")
                data = load(obj.variable_path,obj.variable_name);
                matrix = data.(obj.variable_name);
            else
                matrix = []; %#ok<NASGU>
                error("Cannot locate matrix")
            end
        end
        %----
        function varargout = load_inputs(varargin)
            num_inputs = size(varargin,2);
            varargout = cell(1,num_inputs);
            for iInput = 1:num_inputs   
                input = varargin{iInput};
                if class(input) == "Large_Matrix_Pointer"
                    input = input.load();
                end
                varargout{iInput} = input;
            end

            if nargout == 1
                varargout = {varargout};
            end

        end

        %-------------------
        %-- Overloading
        %----
        function mat_size = size(obj,dim)
            mat_size = obj.matrix_size;
            if nargin == 1
                return
            end
            mat_size = mat_size(dim);
        end
        %----
        function mat = horzcat(mat1,mat2)
            error()
            % mat1_data = mat1.load;
            % mat2_data = mat2.load;
            % mat_data = horzcat(mat1_data,mat2_data);
            % mat = ma1
            % mat1.save(mat_data);
            % mat1.matrix_size = size(mat_data);

        end
        % function sub_matrix = subsindex(obj,index)
        %     matrix = load(obj);
        %     sub_matrix = matrix(index);
        % end
        % %----
        % function output = subsasgn(obj,subscript,value)
        %     matrix = load(obj);
        %     output = subsasgn(matrix,subscript,value);
        % 
        % end
        % %----
        % function output = subsref(obj,subscript)
        %     if subscript.type == "."
        %         output = obj.(subs);
        %         return
        %     end
        %     matrix = load(obj);
        %     output = subsref(matrix,subscript);
        % end

        %----
        function mat_prod = mtimes(a,b)
            [a,b] = load_inputs(a,b);
            mat_prod = a*b;
        end
        %----
        function mat_prod = mrdivide(a,b)
            [a,b] = load_inputs(a,b);
            mat_prod = a/b;
        end
        %----
        function mat_prod = mldivide(a,b)
            [a,b] = load_inputs(a,b);
            mat_prod = a\b;
        end
        %----
        function mat_transpose = ctranspose(a)
            mat = a.load();
            mat_transpose = mat';
        end
        %----
        function varargout = eigs(varargin)
            num_outputs = nargout;
            outputs = load_inputs(varargin{:});
            [V, D, flag] = eigs(outputs{:});
            outputs = {V,D,flag};
            varargout = outputs(1:num_outputs);
        
        end
    end
end