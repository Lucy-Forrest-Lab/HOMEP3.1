# Name: genfsys.py
# Language: python3
# Libraries: sys, os, datetime
# Description: Generates HOMEP file system
# Author: Edoardo Sarti
# Date: Aug 10 2016

import sys, os, datetime
from support import *

# Library function
def generate_filesystem(install_path):
	# Hardcoded variables
	this_name = 'genfsys'
	indent = " "*len(header(this_name))
	version = 3.1
	other_versions_allowed = True

	# Define folder names
	install_path += '/'
	main_dir = 'HOMEP_' + str(version) + '_' + datetime.datetime.now().strftime("%Y_%m_%d") + '/'
	main_path = install_path + main_dir
	rpdb_dir = 'raw_pdbs/'
	cpdb_dir = 'pdbs/'
	lib_dir = {}
	lib_dir['alpha'] = 'alpha/'
	lib_dir['beta'] = 'beta/'

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

	os.mkdir(main_path + rpdb_dir)
	logmsg = header(this_name) + "Directory to store raw pdbs created: {0}\n".format(main_path + rpdb_dir)
	print(logmsg)
	log += logmsg

	os.mkdir(main_path + cpdb_dir)
	logmsg = header(this_name) + "Directory to store curated pdbs created: {0}\n".format(main_path + cpdb_dir)
	print(logmsg)
	log += logmsg


	for ss in 'alpha', 'beta':
		os.mkdir(main_path + lib_dir[ss])
		logmsg = header(this_name) + "Directory to store " + ss + " superfamilies created: {0}\n".format(main_path + lib_dir[ss])
		print(logmsg)
		log += logmsg

	write_log(this_name, log)

	# Compiling output
	locations = {'FSYS' : {}, 'OPT' : {}}
	locations['FSYS']['installpath'] = install_path
	locations['FSYS']['mainpath'] = main_path
	locations['FSYS']['main'] = main_dir
	locations['FSYS']['rpdb'] = rpdb_dir
	locations['FSYS']['cpdb'] = cpdb_dir
	locations['FSYS']['alpha'] = lib_dir['alpha']
	locations['FSYS']['beta'] = lib_dir['beta']

	locations_filename = locations['FSYS']['mainpath'] + '.locations.dat'
	locations_file = open(locations_filename, 'w')
	for x in list(locations['FSYS'].keys()):
		locations_file.write("{0}\t\t{1}\t\t{2}\n".format('FSYS', x, locations['FSYS'][x]))
	locations_file.close()

	return locations


def filesystem_info(main_path):
	this_name = 'chfsys'
	indent = " "*len(header(this_name))
	version = 3.1

	if not os.path.exists(main_path):
		logmsg = header(this_name) + "ERROR: Main directory {0} not found.".format(main_path)
		write_log(this_name, logmsg)
		raise NameError(logmsg)

	locations_filename = main_path + '/' + '.locations.dat'
	if not os.path.exists(locations_filename):
		logmsg = header(this_name) + "ERROR: File {0} not found. Filesystem corrupted.".format(locations_filename)
		write_log(this_name, logmsg)
		raise NameError(logmsg)

	locations = {}
	locations_file = open(locations_filename, 'r')
	text = locations_file.read().split('\n')
	for line in text:
		if line:
			print(line)
			fields = line.split()
			if fields[0] not in locations:
				locations[fields[0]] = {}
			locations[fields[0]][fields[1]] = fields[2]
	return locations
