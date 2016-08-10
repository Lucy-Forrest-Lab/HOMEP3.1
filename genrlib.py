# Name: genrlib.py
# Language: python3
# Libraries:  
# Description: Generates HOMEP raw pdb library
# Author: Edoardo Sarti
# Date: Aug 10 2016

import sys, os, datetime

# Support functions
def write_log(name, text):
	log_filename = name + '.log'
	log_file = open(log_filename, 'w')
	log_file.write(text)
	log_file.close()

def header(name):
	return "[" + name + " " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "] "


# Library function
def generate_raw_pdb_library(names, pdbtm_file_path):
	# Hardcoded variables
	this_name = 'genrlib'
	indent = " "*len(header(this_name))
	version = 3.1

	# Checks
	for path_name in [names['installpath']+x for x in names if x != 'installpath']:
		if not os.path.exists(path_name):
			logmsg = header(this_name) + "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name)
			write_log(this_name, logmsg)	
			raise NameError(logmsg)

	# 
	pdbtm_file = open(pdbtm_file_path, 'r')
	text = pdbtm_file.read().split('\n')
	pdbtm_file.close()
	for line in text:
		fields 
