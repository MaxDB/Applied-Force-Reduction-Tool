function send_email(subject,body)
% sends an email (tested for protonmail)


email_data_path = get_project_path + "\settings\email_data.mat";
script_path = get_project_path + "\code\project_helpers\";
if ~isfile(email_data_path)
    warning("email credentials not found")
    return
end

load(email_data_path,"email_data");

email_data.subject = subject;
email_data.body = body;


send_email_path = "'" + script_path + "send_email_matlab.py'";
pyrunfile(send_email_path,email_data = email_data)
end