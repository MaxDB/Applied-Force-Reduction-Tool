function simulation_update(string_format,values)
subject = sprintf(string_format,values);
body = string(datetime);
send_email(subject,body);
end