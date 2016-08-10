# Name: genfsys.py
# Language: python3
# Libraries: sys, os, datetime
# Description: Generates HOMEP file system
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
def generate_filesystem(install_path):
	# Hardcoded variables
	this_name = 'genfsys'
	indent = " "*len(header(this_name))
	version = 3.1
	other_versions_allowed = True

	# Define folder names
	install_path += '/'
	main_dir = 'HOMEP_' + str(version) + '_' + datetime.datetime.now().strftime("%Y_%m_%d") 
	main_path = install_path + '/' + main_dir
	rpdb_dir = 'raw_pdbs'
	cpdb_dir = 'pdbs'
	lib_dir = {}
	lib_dir['alpha'] = 'alpha'
	lib_dir['beta'] = 'beta'

	# Checks
	if not os.path.exists(install_path):
		logmsg = header(this_name) + "ERROR: The installation directory path {0} does not exist. Please specify and existing path.".format(install_path)
		write_log(this_name, logmsg)	
		raise NameError(logmsg)
	if os.path.exists(main_path):
		logmsg = header(this_name) + "ERROR: In the installation directory path {0} there already is a folder named {1}\n".format(install_path, main_dir)
		write_log(this_name, logmsg)
		raise NameError(logmsg)
	for filename in os.listdir(install_path):
		if filename[0:5] == 'HOMEP' and not other_versions_allowed:
			logmsg = (header(this_name) + "ERROR: In the installation directory path {0} there are other versions of HOMEP.\n".format(install_path) +
			                     indent + "       If you want to continue, you have to set the internal variable other_versions_allowed as True.")
			write_log(this_name, logmsg)
			raise NameError(this_name, logmsg)

	# Generate filesystem
	log = ""

	os.mkdir(main_path)
	logmsg = header(this_name) + "Main directory created: {0}\n".format(main_path)
	print(logmsg)
	log += logmsg

	os.mkdir(main_path + '/' + rpdb_dir)
	logmsg = header(this_name) + "Directory to store raw pdbs created: {0}\n".format(main_path + '/' + rpdb_dir)
	print(logmsg)
	log += logmsg

	os.mkdir(main_path + '/' + cpdb_dir)
	logmsg = header(this_name) + "Directory to store curated pdbs created: {0}\n".format(main_path + '/' + cpdb_dir)
	print(logmsg)
	log += logmsg


	for ss in 'alpha', 'beta':
		os.mkdir(main_path + '/' + lib_dir[ss])
		logmsg = header(this_name) + "Directory to store " + ss + " superfamilies created: {0}\n".format(main_path + '/' + lib_dir[ss])
		print(logmsg)
		log += logmsg

	write_log(this_name, log)

	# Compiling output
	locations = {}
	locations['installpath'] = install_path
	locations['main'] = main_dir
	locations['rpdb'] = rpdb_dir
	locations['cpdb'] = cpdb_dir
	locations['alpha'] = lib_dir['alpha']
	locations['beta'] = lib_dir['beta']

	return locations
