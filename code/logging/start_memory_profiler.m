function start_memory_profiler(varargin)
SAMPLE_DEALY = 1;
log_path = pwd + "\data\logs";
memory_file = "memory.txt";
stop_file = "memory.stop";


%----------------
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

show_console = 1;

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "display"
            show_console = keyword_values{arg_counter};
        case {"colour","color"}
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-------------------------------------------------------------------------%


if ~isfolder(log_path)
    mkdir(log_path)
end

%check for file
memory_path = log_path + "\" + memory_file;
stop_path = log_path + "\" + stop_file;
if isfile(memory_path) && ~isfile(stop_path)
    stop_memory_profiler
    pause(SAMPLE_DEALY*1.1)
end

if isfile(memory_path)
    delete(memory_path)
end

if isfile(stop_path)
    delete(stop_path)
end


%run profiler
memory_log_path = convertStringsToChars(log_path);
script_path =  convertStringsToChars(get_project_path + "\code\logging");

if ~ischar(show_console)
    show_console = convertStringsToChars(string(show_console));
end

start_script_path = [script_path,'\memory_profiler.ps1'];
powershell_command = ['Start-Process powershell \"-File \"\"',start_script_path,'\"\" -\"\"log_path\"\" \"\"',memory_log_path,'\"\" -\"\"display\"\" ',show_console,' \"'];
cmd_command = ['powershell -Command "& {',powershell_command,'}"'];
system(cmd_command);

end