function abaqus_detected = abaqus_check()

command = "abaqus -whereami";
[status,msg] = system(command);

abaqus_detected = ~status;

msg_lines = splitlines(msg);
abaqus_version = msg_lines{1};

if abaqus_detected
    cmd_msg = fprintf('%s detected\n',abaqus_version);
else
    cmd_msg_1 = fprintf(2,"Warning: Cannot detect Abaqus: ");
    cmd_msg_2 = fprintf('Please ensure Abaqus is on the system path and accessible with the alias "abaqus"\n');
    cmd_msg = cmd_msg_1 + cmd_msg_2;
end

end_line = repelem('_',cmd_msg - 1);
disp(end_line)
end