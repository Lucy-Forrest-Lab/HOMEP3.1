# Name: support.py
# Language: python3
# Libraries:
# Description: Support functions used by modules
# Author: Edoardo Sarti
# Date: Aug 10 2016

import datetime

# Support functions
def write_log(name, text):
	log_filename = name + '.log'
	log_file = open(log_filename, 'w')
	log_file.write(text)
	log_file.close()

def header(name):
	return "[" + name + " " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "] "

