console_size = matlab.desktop.commandwindow.size;
row_length = console_size(1)-1;
bar_line = repelem('_',console_size(1)-1);



afr_tool_parts = {...
'    ___    __________     ______            __',...
'   /   |  / ____/ __ \   /_  __/___  ____  / /',...
'  / /| | / /_  / /_/ /    / / / __ \/ __ \/ / ',...
' / ___ |/ __/ / _, _/    / / / /_/ / /_/ / /  ',...
'/_/  |_/_/   /_/ |_|    /_/  \____/\____/_/   '...
};
%generated using https://patorjk.com/software/taag/ with the 'Slant' font



name_length = length(afr_tool_parts{1});
max_padding = row_length - name_length;
padding_length = floor(max_padding/2);
padding = repelem(' ',padding_length);
afr_tool_parts = cellfun(@(line) [padding,line,padding] ,afr_tool_parts,"UniformOutput",false);

                                     
afr_tool = sprintf('%s\n',afr_tool_parts{:});
afr_tool(end) = [];

fprintf('\n')
disp(bar_line)
disp(afr_tool)
disp(bar_line)


%-----------------
current_version_text = "Installed version: " + get_version();
disp(current_version_text)
line_length = fprintf("Checking for updates...");

is_new_version = check_version;
fprintf(repmat('\b',1,line_length))

if isempty(is_new_version)
    new_version_text = "Unable to check current version";
else
    if is_new_version
        new_version_text = "New version available at https://github.com/MaxDB/Applied-Force-Reduction-Tool"; 
    else
        new_version_text = "Version up to date";
    end
end
disp(new_version_text)

version_length = max(strlength(new_version_text),strlength(current_version_text));
end_line = repelem('_',version_length);
disp(end_line)

%-------------------
function version = get_version
    project_path = get_project_path;
    version_path = project_path + "\settings\version_number.txt";
    if ~isfile(version_path)
        version = "Error - could not find version number";
        return
    end

    version_id = fopen(version_path);
    try
        version = textscan(version_id,"%s");
        version = convertCharsToStrings(version{1}{1});
    catch exception
        fclose(version_id);
        rethrow(exception)
    end
    fclose(version_id);
end


function is_new_version=check_version()
version_path = "settings\version_number.txt";
try
    [exit_code(1),cmd_output] = system("git fetch origin main"); %#ok<*ASGLU>
    [exit_code(2),cmd_output] = system("git diff origin/main " + version_path);
    if any(exit_code ~= 0) 
        is_new_version = ~isempty(cmd_output);
    else
        is_new_version = logical([]);
    end
catch
    is_new_version = logical([]);
end


end