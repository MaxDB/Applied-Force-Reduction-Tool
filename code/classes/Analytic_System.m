classdef Analytic_System
    %System with known equations of motion
    properties
        system_name 
        linear_mass
        linear_damping
        linear_stiffness
        nonlinear_restoring_force
        nonlinear_restoring_force_jacobian
        external_force
        potential_energy

        eigenvectors
        eigenvalues

        modal_restoring_force
        modal_stiffness
        modal_potential_energy
    end
    methods
        function obj = Analytic_System(name,eom,varargin)
            
            %Optional argumanents
            num_args = length(varargin);
            if mod(num_args,2) == 1
                error("Invalid keyword/argument pairs")
            end
            keyword_args = varargin(1:2:num_args);
            keyword_values = varargin(2:2:num_args);
            
            save_system = false; %save system locally

            for arg_counter = 1:num_args/2
                switch keyword_args{arg_counter}
                    case "save"
                        save_system = true;
                        save_path = keyword_values{arg_counter};
                        continue
                    otherwise
                        error("Invalid keyword: " + keyword_args{arg_counter})
                end
            end
            %%%%%

            obj.system_name = name;

            if isfield(eom,'M')
                obj.linear_mass = eom.M;
            end
            if isfield(eom,'K')
                obj.linear_stiffness = eom.K;
            end
            if isfield(eom,'f')
                obj.nonlinear_restoring_force = eom.f;
                obj.nonlinear_restoring_force_jacobian = obj.get_jacobian(eom.f);
            end
            if isfield(eom,'F')
                obj.external_force = eom.F;
            end
            if isfield(eom,'C')
                obj.linear_damping = eom.C;
            end
            if isfield(eom,'V')
                obj.potential_energy = eom.V;
            end

            [evecs,evals] = eig(obj.linear_stiffness,obj.linear_mass,"vector");
            obj.eigenvectors = evecs;
            obj.eigenvalues = evals;

            [fq,dfqdq] = obj.get_modal_restoring_force;
            obj.modal_restoring_force = fq;
            obj.modal_stiffness = dfqdq;
            obj.modal_potential_energy = obj.get_modal_potential_energy;

            if save_system
                obj.save_class_instance(save_path)
            end
        end
        %-----------------------------------------------------------------%
        function save_class_instance(analytic_eom,save_path)
            if isfolder(save_path)
                rmdir(save_path,"s")
            end
            mkdir(save_path)
            save(save_path + "\" + analytic_eom.system_name,"analytic_eom")
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function static_equation = get_static_equation(obj)
            static_equation = @(x) obj.linear_stiffness*x + obj.nonlinear_restoring_force(x);
        end
        %-----------------------------------------------------------------%
        function potential_equation = get_potential_equation(obj)
            if ~isempty(obj.potential_energy)
                potential_equation = obj.potential_energy;
                return
            end

            % Integrate work
            % Parameterise in terms of t first
            static_equation = obj.get_static_equation;
        end
        %-----------------------------------------------------------------%
        function modal_potential = get_modal_potential_energy(obj)
            modal_potential = @(q) obj.potential_energy(obj.eigenvectors*q);
        end
        %-----------------------------------------------------------------%
        function J = get_jacobian(obj,f)
            num_dofs = size(obj.linear_mass,1);
            x = sym("x",[num_dofs,1]);
            J_sym = jacobian(f(x),x);
            J = matlabFunction(J_sym);
        end
        %-----------------------------------------------------------------%
        function [fq,dfqdq] = get_modal_restoring_force(obj)
            evecs = obj.eigenvectors;
            fx = obj.nonlinear_restoring_force;
            fq = @(q) diag(obj.eigenvalues)*q + evecs'*fx(evecs*q); 

            num_dofs = size(evecs,1);
            q = sym("q",[num_dofs,1]);
            J_sym = jacobian(fq(q),q);
            dfqdq = matlabFunction(J_sym);
        end
        %-----------------------------------------------------------------%
        function Eom_Input = get_solver_inputs(obj,type)
            switch type
                case "free"
                    Eom_Input.modal_restoring_force = obj.modal_restoring_force;
                    Eom_Input.modal_stiffness = obj.modal_stiffness;
                    Eom_Input.modal_potential = obj.modal_potential_energy;
                case "forced"

            end
            
        end
    end
end