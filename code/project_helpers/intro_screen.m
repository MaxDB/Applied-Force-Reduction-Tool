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

disp(bar_line)
disp(afr_tool)
disp(bar_line)


%-----------------
current_version_text = "Installed version: " + get_version();
is_new_version = check_version;
disp(current_version_text)
disp(available_version_text)



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
