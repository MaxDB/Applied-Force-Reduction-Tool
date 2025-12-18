function simulation_update(string_format,values)
subject = sprintf(string_format,values);
body = string(datetime);

try
    send_email(subject,body);
catch matlab_error
    warning(matlab_error.message)
    warning("Email not sent")
end