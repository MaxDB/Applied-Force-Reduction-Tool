function logger(log_message,message_level)
timestamp_format = "HH:mm:ss.SS";
timestamp = datetime;
timestamp.Format = timestamp_format;


worker = getCurrentWorker();


log_path = "data\logs\";
log_level = "data\log_level.mat";
log_file = log_path + "log.txt";
if ~isempty(worker)
    worker_id = worker.ProcessId;
    lock_file = log_path + "log_" + worker_id + ".lck";
else
    lock_file = log_path + "log.lck";
end



load(log_level,"logging_level")

switch message_level
    case 1
        prefix = "";
    case 2
        prefix = "\t";
    case 3
        prefix = "\t\t";
    case 4
        prefix = "\t\t\t";
    otherwise
        prefix = "";
end

suffix = "\n";
log_message = prefix + log_message + suffix;


if message_level == 1
    log_break = "---------------------------------------------------------";
    logger(log_break,2)
end

if message_level <= logging_level
    fprintf(log_message)
end

log_message = string(timestamp) + ">>\t" + log_message;

lock_counter = 0;
while ~isempty(dir("*.lck"))
    pause(0.1);
    lock_counter = lock_counter + 1;
    if lock_counter > 10
        error("Log locked")
    end
end


lock_id = fopen(lock_file,'w');
fclose(lock_id);
log_id = fopen(log_file,"a");
try
    fprintf(log_id,log_message);
    fclose(log_id);
catch
    warning("Could not write to log")
    fclose(log_id);
end
delete(lock_file)

if message_level == 1
    logger("",2)
end