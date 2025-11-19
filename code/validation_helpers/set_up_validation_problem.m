function [h_terms,reduced_eom,h_solver,Validation_Analysis_Inputs] = set_up_validation_problem(Validation_Rom,Validation_Opts,Solution,Validated_BB_Settings)

Solution_Type = Solution.Solution_Type;
orbit_type = Solution_Type.orbit_type;

switch orbit_type
    case "free"
        Eom_Input = Validation_Rom.get_solver_inputs("coco_backbone");
        reduced_eom = @(t,z,zeta) coco_eom(t,z,zeta,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data);
        %
        Validation_Input = Validation_Rom.get_solver_inputs("h_prediction");
        h_terms = @(r,r_dot,r_ddot) get_h_error_terms(r,r_dot,r_ddot,Validation_Input);

        if nargout == 1
            return
        end

        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("h_analysis","additional_output",Validated_BB_Settings.Additional_Output);
    case "forced"
        
        Nonconservative_Input = Solution.get_nonconservative_input(Validation_Rom.Model);
        amp = Nonconservative_Input.amplitude;
        Eom_Input = Validation_Rom.get_solver_inputs("coco_frf",Nonconservative_Input);
        reduced_eom = @(t,z,T) coco_forced_eom(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);
        
        Validation_Input = Validation_Rom.get_solver_inputs("forced_h_prediction",Nonconservative_Input);
        h_terms = @(t,r,r_dot,r_ddot,period) get_forced_h_error_terms(t,r,r_dot,r_ddot,amp,period,Validation_Input);

        if nargout == 1
            return
        end

        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("forced_h_analysis",Nonconservative_Input);
end


switch Validation_Opts.validation_algorithm
    case "h_time"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_time_solution(h_terms,t0,omega,num_harmonics);
    case "h_frequency"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_harmonic_balance(h_terms,t0,omega,num_harmonics);
    case "h_infinite_determinant"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_infinite_determinant(h_terms,t0,omega,num_harmonics);
end