function start_crash_monitoring()
% sends an email if matlab crashes (protonmail only)


email_data_path = get_project_path + "\settings\email_data.mat";
script_path = get_project_path + "\code\project_helpers\";
if ~isfile(email_data_path)
    warning("email credentials not found")
    return
end

load(email_data_path,"email_data");


if double(version(1:2)) >= 25
    matlab_pid = matlabProcessID;
else
    matlab_pid = feature('getpid');
end

monitor_matlab_path = "'" + script_path + "monitor_matlab.py'";
pyrunfile(monitor_matlab_path,matlab_pid = matlab_pid,email_data = email_data)
% email_data.body = "matlab test";
% email_data.subject = "matlab test 1";
% 
% send_email_path = "'" + script_path + "send_email_matlab.py'";
% pyrunfile(send_email_path,email_data = email_data)

end
