# Name: support.py
# Language: python3
# Libraries: datetime
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


def raise_error(name, text):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(this_name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	write_log(name, logmsg)
	raise NameError(logmsg)


def print_log(name, text):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(this_name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	print(logmsg)
	return logmsg
