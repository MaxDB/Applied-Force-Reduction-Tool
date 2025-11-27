from datetime import datetime
import time
import psutil

import sys
sys.path.append("C:\\Users\\gg19546\\Documents\\PhD\\Year 2\\Software\\Applied Force Reduction\\code\\project_helpers")
from send_email import send_email

#expects matlab_pid and email_data in input

sender = email_data['sender']
recipient = email_data['reciever']
port_number = int(email_data['port'])
password = email_data['password']

check_rate = 1

while True:
    time.sleep(check_rate)
    if not psutil.pid_exists(matlab_pid):
        subject = "Simulation crashed: " + str(datetime.now())
        body = ":("
        send_email(sender,recipient,subject,body,port_number,password)
        break
    print("test")