function [num_sep_lines,additional_job_input_lines] = estimate_input_lines(num_dofs,Static_Opts)
num_loadcases = Static_Opts.num_loadcases;
loadcase_lines = num_loadcases*num_dofs;

switch Static_Opts.additional_data
    case "perturbation"
        validation_lines = Static_Opts.num_validation_modes*loadcase_lines;
    otherwise
        validation_lines = 0;
end


%ignore sep reset for now
num_sep_lines = validation_lines + loadcase_lines;

additional_job_input_lines = num_dofs; %(divide by dimension
end