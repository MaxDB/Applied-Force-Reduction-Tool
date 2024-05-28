function logger(log_message,message_level)
load("data\log_level.mat","logging_level")

switch message_level
    case 1
        prefix = "";
    case 2
        prefix = "\t";
    case 3
        prefix = "\t\t";
    otherwise
        prefix = "";
end

suffix = "\n";
log_message = prefix + log_message + suffix;


if message_level == 1
    log_break = "---------------------------------------------------------";
    logger(log_break,2)
end

log_id = fopen("data\logs\log.txt","a");
fprintf(log_id,log_message);
fclose(log_id);

if message_level <= logging_level
    fprintf(log_message)
end

if message_level == 1
    logger("",2)
end