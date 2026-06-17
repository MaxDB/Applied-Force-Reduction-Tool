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
        
        Nonconservative_Input = Forced_Solution.get_nonconservative_input(Solution.Force_Data,Solution.Damping_Data,Validation_Rom);
        amp = Nonconservative_Input.amplitude;
        Eom_Input = Validation_Rom.get_solver_inputs("coco_frf","additional_input",Nonconservative_Input);
        reduced_eom = @(t,z,T) coco_forced_eom(t,z,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);
        
        Validation_Input = Validation_Rom.get_solver_inputs("forced_h_prediction","additional_input",Nonconservative_Input);
        h_terms = @(t,r,r_dot,r_ddot,period) get_forced_h_error_terms(t,r,r_dot,r_ddot,amp,period,Validation_Input,1);

        if nargout == 1
            return
        end

        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("forced_h_analysis","additional_input",Nonconservative_Input,"additional_output",Validated_BB_Settings.Additional_Output);
    case "frf_to_bb"
        Nonconservative_Input = Solution.get_nonconservative_input(Validation_Rom);
        amp = Nonconservative_Input.amplitude;
        frequency = Nonconservative_Input.frequency;
        T = 2*pi/frequency;
        Eom_Input = Validation_Rom.get_solver_inputs("coco_frf","additional_input",Nonconservative_Input);
        reduced_eom = @(t,z,epsilon) coco_frf2bb_eom(t,z,epsilon,amp,T,Eom_Input.input_order,Eom_Input.Force_Data,Eom_Input.Disp_Data,Eom_Input.Damping_Data,Eom_Input.Applied_Force_Data);

        Validation_Input = Validation_Rom.get_solver_inputs("forced_h_prediction","additional_input",Nonconservative_Input);
        h_terms = @(t,r,r_dot,r_ddot,epsilon) get_forced_h_error_terms(t,r,r_dot,r_ddot,amp,T,Validation_Input,epsilon);

        if nargout == 1
            return
        end
       
        Validation_Analysis_Inputs = Validation_Rom.get_solver_inputs("forced_h_analysis","additional_input",Nonconservative_Input,"additional_output",Validated_BB_Settings.Additional_Output);
end


switch Validation_Opts.validation_algorithm
    case "h_time"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_time_solution(h_terms,t0,omega,num_harmonics);
    case "h_frequency"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_harmonic_balance(h_terms,t0,omega,num_harmonics);
    case "h_infinite_determinant"
        h_solver = @(h_terms,t0,omega,num_harmonics) h_infinite_determinant(h_terms,t0,omega,num_harmonics);
end