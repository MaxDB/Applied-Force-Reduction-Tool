function [memory_data,sample_period,date] = get_free_memory
log_path = pwd + "\data\logs";
memory_file = "memory.txt";

memory_path = log_path + "\" + memory_file;
if ~isfile(memory_path)
    error("Cannot find " + memory_path)
end


memory_id = fopen(memory_path);
try
    date = textscan(memory_id,"%D %D",1);
    sample_period = textscan(memory_id,"%f",1);
    memory_data = textscan(memory_id,"%f");
catch exception
    fclose(memory_id);
    rethrow(exception)
end
fclose(memory_id);

memory_data = memory_data{1};
sample_period = sample_period{1};
date = date{1} + timeofday(date{2});
end