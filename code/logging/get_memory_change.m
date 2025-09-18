function memory_change = get_memory_change
log_path = pwd + "\data\logs";
memory_file = "memory.txt";

memory_path = log_path + "\" + memory_file;
if ~isfile(memory_path)
    error("Cannot find " + memory_path)
end

memory_id = fopen(memory_path);
try
    memory_data = textscan(memory_id,"%f");
catch exception
    fclose(memory_id);
    rethrow(exception)
end
fclose(memory_id);

memory_data = memory_data{1};
memory_change = max(memory_data) - min(memory_data);
end