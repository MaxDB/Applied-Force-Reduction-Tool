function stop_memory_profiler
log_path = pwd + "\data\logs";
stop_file = "memory.stop";
memory_file = "memory.txt";


stop_path = log_path + "\" + stop_file;
memory_path = log_path + "\" + memory_file;

if ~isfile(memory_path)
    return
end

cmd_command = "type nul >""" + stop_path + """";
system(cmd_command)

end