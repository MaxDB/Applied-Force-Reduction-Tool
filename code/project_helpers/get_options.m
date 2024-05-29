function options = get_options(options_name)
project_path = get_project_path;
options_path = project_path + "\settings\" + options_name + "_options.txt";
opts_id = fopen(options_path);
options_data = textscan(opts_id,'%s','delimiter','\n');
fclose(opts_id);
options_data = options_data{1,1};

num_lines = size(options_data,1);
options = cell(num_lines,3);
option_counter = 0;
for iLine = 1:num_lines
    line = strip(options_data{iLine,1});
    if startsWith(line,"#")
        continue
    end
    if isempty(line)
        continue
    end

    line_data = split(line,",");
    option_name = strip(line_data{1,1});
    option_value = strip(line_data{2,1});
    if startsWith(option_value,'"')
        option_value = string(strip(option_value,"both",'"'));
    elseif startsWith(option_value,'[')
        option_value = string(strip(option_value,"left",'['));
        option_value = string(strip(option_value,"right",']'));
        option_value = split(option_value,";");
        option_value = str2double(option_value);
    else
        option_value = str2double(option_value);
    end
    option_counter = option_counter + 1;
    option_description = strip(line_data{3,1});
    options(option_counter,:) = {option_name,option_value,option_description};
end

options = options(1:option_counter,:);
end