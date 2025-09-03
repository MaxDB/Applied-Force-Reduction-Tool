function logger(log_message,message_level)
timestamp_format = "HH:mm:ss.SS";
timestamp = datetime;
timestamp.Format = timestamp_format;

load("data\log_level.mat","logging_level")

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

log_id = fopen("data\logs\log.txt","a");
fprintf(log_id,log_message);
fclose(log_id);

if message_level == 1
    logger("",2)
end