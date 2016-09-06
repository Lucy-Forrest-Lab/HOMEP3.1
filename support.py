# Name: support.py
# Language: python3
# Libraries: datetime
# Description: Support functions used by modules
# Author: Edoardo Sarti
# Date: Aug 10 2016

import datetime
import shutil
import os


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
	logmsg = header(name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	write_log(name, logmsg)
	raise NameError(logmsg)


def print_log(name, text):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	print(logmsg)
	return logmsg


def archive_old_file(locations, filenames):
	if type(filenames) == str:
		filenames = [filenames]
	for filename in filenames:
		filename_noext = os.path.basename(os.path.splitext(filename)[0])
		extension = os.path.splitext(filename)[1]
		one_hundred_files = False
		for i in range (0,100):
			new_filename = locations['FSYS']['mainpath'] + '.old/' + filename_noext + '_' + str(i) + extension
			if not os.path.exists(new_filename):
				shutil.move(filename, new_filename)
				break
			if i == 99:
				one_hundred_files = True
		if one_hundred_files:
			for i in range(1, 100):
				old_filename = locations['FSYS']['mainpath'] + '.old/' + filename_noext + '_' + str(i) + extension
				new_filename = locations['FSYS']['mainpath'] + '.old/' + filename_noext + '_' + str(i-1) + extension
				shutil.move(old_filename, new_filename)
			new_filename = locations['FSYS']['mainpath'] + '.old/' + filename_noext + '_' + str(99) + extension
			shutil.move(filename, new_filename)
