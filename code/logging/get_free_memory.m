function [memory_data,memory_duration,approx_sample_period] = get_free_memory
log_path = pwd + "\data\logs";
memory_file = "memory.txt";
SEPERATOR = "--";

memory_path = log_path + "\" + memory_file;
if ~isfile(memory_path)
    error("Cannot find " + memory_path)
end


memory_id = fopen(memory_path);
date_format = "%{dd/MM/uuuu}D";
time_format = "%{HH:mm:ss}D";
try
    start_date = textscan(memory_id,date_format + " " + time_format,1);
    approx_sample_period = textscan(memory_id,"%f",1);
    textscan(memory_id,SEPERATOR,1);
    memory_data = textscan(memory_id,"%f");
    textscan(memory_id,SEPERATOR,1);
    end_date = textscan(memory_id,date_format + " " + time_format,1);
catch exception
    fclose(memory_id);
    rethrow(exception)
end
fclose(memory_id);

memory_data = memory_data{1};
approx_sample_period = approx_sample_period{1};

start_date = start_date{1} + timeofday(start_date{2});
end_date = end_date{1} + timeofday(end_date{2});

memory_duration = seconds(end_date-start_date);
end