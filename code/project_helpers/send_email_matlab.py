sender = email_data['sender']
recipient = email_data['reciever']
subject = email_data['subject']
body = email_data['body']   
port_number = int(email_data['port'])
password = email_data['password']


import smtplib 
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

def send_email(sender,recipient,subject,body,port_number,password):
    msg = MIMEMultipart()
    msg['From'] = sender
    msg['To'] = recipient
    msg['Subject'] = subject
    message = body
    msg.attach(MIMEText(message))
    mailserver = smtplib.SMTP('localhost',port_number)
    mailserver.login(sender, password)
    mailserver.sendmail(sender,recipient,msg.as_string())
    mailserver.quit()


send_email(sender,recipient,subject,body,port_number,password)