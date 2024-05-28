classdef Sparse_Stiffness_Array
    %Stores an array of sparse matrices with the same sparsity pattern
    properties
        nonzero_indicies
        nonzero_data

        dimensions
        num_load_cases

        symmetrical
    end
    methods
        function obj = Sparse_Stiffness_Array(dimensions,sparsity_pattern,symmetry,num_loadcases)
            obj.dimensions = dimensions;
            obj.symmetrical = symmetry;
            obj.nonzero_indicies = sparsity_pattern;
            obj.nonzero_data = zeros(size(sparsity_pattern,1), dimensions(3));
            if nargin == 3
                obj.num_load_cases = 0;
            else
                obj.num_load_cases = num_loadcases;
            end
        end
        %-----------------------------------------------------------------%
        function obj = add_load_case(obj,K)
            obj.num_load_cases = obj.num_load_cases + 1;
            obj.nonzero_data(:,obj.num_load_cases) = K(:,end);
        end
        %-----------------------------------------------------------------%
        function K = get_matrix(obj,load_index)
            nz_indicies = obj.nonzero_indicies;
            nz_data = obj.nonzero_data;
            nz_data_point = nz_data(:,load_index);
            K = sparse(nz_indicies(:,1),nz_indicies(:,2),nz_data_point);
        end
        %-----------------------------------------------------------------%
        function K_i = get_component(obj,varargin)
            if nargin == 2
                %linear index
                lin_index = varargin{1,1};
                K_i = obj.nonzero_data(lin_index,:);
            elseif nargin == 3
                %coordinate
                
            end
        end
        %-----------------------------------------------------------------%
        function nz_data = match_sparsity_pattern(obj,K)
            nz_indicies = obj.nonzero_indicies;
            num_rows = size(nz_indicies,1);
            nz_data = zeros(num_rows,1);
            for iRow = 1:num_rows
                nz_data(iRow,1) = K(nz_indicies(iRow,1),nz_indicies(iRow,2));
            end
           
        end
        %-----------------------------------------------------------------%
        
        %-- Overloading
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function sz = size(obj,dim)
            sz = obj.dimensions;
            if nargin == 2
                sz = sz(dim);
            end
        end
        %-----------------------------------------------------------------%
        function obj_1 = horzcat(obj_1,obj_2)
            if ~isequal(obj_1.nonzero_indicies,obj_2.nonzero_indicies)
                error("Cannot concatenate: different sparsity patterns")
            end
            if ~isequal(obj_1.dimensions(1:(end-1)),obj_2.dimensions(1:(end-1)))
                error("Cannot concatenate: different sizes")
            end
            
            num_loads = obj_1.num_load_cases + obj_2.num_load_cases;
            obj_1.num_load_cases = num_loads;
            obj_1.dimensions(end) = num_loads;
            obj_1.nonzero_data = [obj_1.nonzero_data,obj_2.nonzero_data];
        end
        %-----------------------------------------------------------------%
        function obj_1 = cat(dim,obj_1,obj_2)
            if ~isequal(obj_1.nonzero_indicies,obj_2.nonzero_indicies)
                error("Cannot concatenate: different sparsity patterns")
            end
            if ~isequal(obj_1.dimensions(1:(end-1)),obj_2.dimensions(1:(end-1)))
                error("Cannot concatenate: different sizes")
            end
            
            num_loads = obj_1.num_load_cases + obj_2.num_load_cases;
            obj_1.num_load_cases = num_loads;
            obj_1.dimensions(end) = num_loads;
            obj_1.nonzero_data = [obj_1.nonzero_data,obj_2.nonzero_data];
        end
        %-----------------------------------------------------------------%
        function obj = subsasgn(obj,Assignment,stiffness)
            type = Assignment.type;
            subs = Assignment.subs;
            switch type
                case "{}"
                    obj.nonzero_data(:,subs{1,1}) = stiffness(:,end);
            end
        end
    end

end