function log_data = read_log(log_lines)
log_path = "data\logs\log.txt";
timestamp_format = "HH:mm:ss.SS";

if isfile(log_path)
    log_id = fopen(log_path);
else
    error("log doesn't exist")
end


try
    log_file = textscan(log_id,"%s %s","Delimiter",{'>>','\n'});
catch
    fclose(log_id);
    error("Cannot read log")
end
fclose(log_id);


%------------------
log_times = log_file{1};
log_text = log_file{2};


num_outputs = size(log_lines,1);
log_timespan = zeros(num_outputs,1);
for iOutput = 1:num_outputs
    line_start = log_lines(iOutput,1);
    log_index = find(startsWith(log_text,line_start),1);
    if ~isempty(log_index)
        timestamp_start = log_times{log_index};
        start_time = datetime(timestamp_start,"InputFormat",timestamp_format,"Format",timestamp_format);
        log_text(1:log_index) = [];
        log_times(1:log_index) = [];
    else
        start_time = nan;
    end

    line_end = log_lines(iOutput,2);
    log_index = find(startsWith(log_text,line_end),1);
    if ~isempty(log_index)
        timestamp_end = log_times{log_index};
        end_time = datetime(timestamp_end,"InputFormat",timestamp_format,"Format",timestamp_format);
    else
        end_time = nan;
    end
    
    if class(start_time) == "datetime" && class(end_time) == "datetime"
        time_diff = seconds(end_time - start_time);
        if time_diff < 0
            time_diff = time_diff + 24*60*60;
        end
    else
        time_diff = nan;
    end
    log_timespan(iOutput) = time_diff;
    
end



log_data = log_timespan;
end


