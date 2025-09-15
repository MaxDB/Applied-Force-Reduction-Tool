function logger(log_message,message_level,varargin)
timestamp_format = "HH:mm:ss.SS";

%-------------------------------------------------------------------------%
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

print_timestamp = 1;
for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "timestamp"
            print_timestamp = keyword_values{arg_counter};
        otherwise
            error("Invalid keyword: " + keyword_args{arg_counter})
    end
end
%-----



worker = getCurrentWorker();


log_path = "data\logs\";
log_level = "data\log_level.mat";
log_file = log_path + "log.txt";
if ~isempty(worker)
    worker_id = worker.ProcessId;
    lock_file = log_path + "log_" + worker_id + ".lck";
else
    worker_id = 0;
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



lock_files_path = "data\logs\*.lck";
lock_counter = 0;
lock_files = dir(lock_files_path);
while ~isempty(lock_files)
    pause(0.01);
    lock_counter = lock_counter + 1;
    if lock_counter > 100
        error("Log locked")
    end
    lock_files = dir(lock_files_path);
end
lock_id = fopen(lock_file,'w');
fclose(lock_id);

while ~isfile(lock_file)
    pause(0.01);
end

lock_files = dir(lock_files_path);
if length(lock_files) > 1
    min_id = 0;
    lock_counter = 0;
    while worker_id ~= min_id
        lock_file_names = arrayfun(@(x) convertCharsToStrings(x.name),lock_files);
        lock_file_ids = double(extract(lock_file_names,digitsPattern));
        min_id = min(lock_file_ids);
        pause(0.01);
        lock_counter = lock_counter + 1;
        if lock_counter > 100
            error("Log locked")
        end
        lock_files = dir(lock_files_path);
    end
end






log_id = fopen(log_file,"a");
try
    if print_timestamp
        timestamp = datetime;
        timestamp.Format = timestamp_format;
        log_message = string(timestamp) + ">>\t" + log_message;
    end
    if message_level <= logging_level
        fprintf(log_message)
    end
    fprintf(log_id,log_message);
    fclose(log_id);
catch exception_message
    fclose(log_id);
    rethrow(exception_message)
end

delete(lock_file)

while isfile(lock_file)
    pause(0.01);
end

if message_level == 1
    logger("",2)
end