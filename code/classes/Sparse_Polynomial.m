classdef Sparse_Polynomial
    %Stores an array of sparse matrices modelled by polynomials
    properties
        nonzero_indicies
        Nonzero_Polynomials

        dimensions
    end

    methods
        function obj = Sparse_Polynomial(Stiffness_Array,Polynomial_Array)
            obj.dimensions = Stiffness_Array.dimensions(1:(end-1));
            obj.nonzero_indicies = Stiffness_Array.nonzero_indicies;
            obj.Nonzero_Polynomials = Polynomial_Array;
        end
        %-----------------------------------------------------------------%
        function plot_polynomial(obj,ax,plotted_outputs)
            obj.Nonzero_Polynomials.plot_polynomial(ax,plotted_outputs);
        end
        %-----------------------------------------------------------------%
        function output_data = evaluate_polynomial(obj,input_data,output)
            if nargin == 2
                output = 1:size(obj.nonzero_indicies,1);
            end
            output_data = obj.Nonzero_Polynomials.evaluate_polynomial(input_data,output);
        end

        %-- Overloading
        %-----------------------------------------------------------------%
        %-----------------------------------------------------------------%
        function sz = size(obj)
            sz = obj.dimensions;
        end
        %-----------------------------------------------------------------%
    end
end