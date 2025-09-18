function [status,cmd_out] = run_abaqus_job(job,varargin)
%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

num_cpus = 1;
interactive = "off"; %"off" / "on" / level
temp = 1;
for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "num_cpus"
            num_cpus = keyword_values{arg_counter};
        case "interactive"
            interactive = keyword_values{arg_counter};
        case "temp"
            temp = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-----
if isstring(interactive)
    if interactive == "on"
        interactive = 0;
    elseif interactive == "off"
        interactive = inf;
    end
end
%-------------------------------------------------------------------------%
log_level = "data\log_level.mat";
load(log_level,"logging_level")

command = "abaqus job=" + job + " cpus=" + num_cpus + " interactive";

num_workers = get_current_parallel_jobs;

is_parallel = num_workers > 1;
log_start_message = "|---------------------- Abaqus output start ----------------------|";
log_end_message = "|----------------------- Abaqus output end -----------------------|";


if ~is_parallel
    logger(log_start_message,interactive)
end
if temp
    project_directory = pwd;
    cd temp
end
if interactive <= logging_level && ~is_parallel
    [status,cmd_out] = system(command,"-echo");
else
    [status,cmd_out] = system(command);
    while ~isfile(job + ".dat")
        pause(0.1)
    end
end

while isfile(job + ".lck")
    pause(0.1)
end

if temp
    cd(project_directory)
end

if is_parallel
    logger(log_start_message,interactive,"timestamp",0)
    logger(cmd_out,interactive,"timestamp",0)
    logger(log_end_message,interactive,"timestamp",0)
else
    logger(cmd_out,inf,"timestamp",0)
    logger(log_end_message,interactive)
end
end
