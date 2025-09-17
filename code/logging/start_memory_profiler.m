function start_memory_profiler
SAMPLE_DEALY = 1;
log_path = pwd + "\data\logs";
memory_file = "memory.txt";
stop_file = "memory.stop";



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

start_script_path = [script_path,'\memory_profiler.ps1'];
powershell_command = ['Start-Process powershell \"-File \"\"',start_script_path,'\"\" -\"\"log_path\"\" \"\"',memory_log_path,'\"\"\"'];
cmd_command = ['powershell -Command "& {',powershell_command,'}"'];
system(cmd_command);

%make sure to add killing the memoey profiler to project shut down
end