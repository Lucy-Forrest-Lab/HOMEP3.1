# Name: genfsys.py
# Language: python3
# Libraries: sys, os, datetime
# Description: Generates HOMEP file system
# Author: Edoardo Sarti
# Date: Aug 10 2016

import sys, os, datetime, argparse
from support import *

def main_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--install_dir', nargs=1)
	parser.add_argument('-pdbtm', '--pdbtm_file_path', nargs=1)
	parser.add_argument('-s', '--straln_path', nargs=1)
	parser.add_argument('-np', '--number_of_procs', nargs=1)
	parser.add_argument('-ot', '--object_thr', nargs=1)
	parser.add_argument('-ct', '--cluster_thr', nargs=1)
	parser.add_argument('-rf', '--resolution_filter', nargs=1)
	parser.add_argument('-with_nmr', action='store_true')
	parser.add_argument('-with_theoretical', action='store_true')
	parser.add_argument('-ht', '--hole_thr', nargs='?')
	parser.add_argument('-oh', '--output_homep', nargs='?')
	parser.set_defaults(hole_thr = '100')
	parser.set_defaults(output_tab = 'structure_alignments.dat')
	parser.set_defaults(output_homep = 'HOMEP3.1.dat')
	parsed = parser.parse_args()

	options = {}
	for x in sorted(parsed.__dict__):
		if type(parsed.__dict__[x]) == list:
			print(x, parsed.__dict__[x][0])
			options[x] = parsed.__dict__[x][0]
		else:
			print(x, parsed.__dict__[x])
			options[x] = parsed.__dict__[x]

	filters = {'resolution' : float(parsed.resolution_filter[0]),
	           'NMR' : parsed.with_nmr,
	           'THM' : parsed.with_theoretical,
	           'hole_thr' : int(parsed.hole_thr[0])}

	return options, filters



# Library function
def generate_filesystem():
	# Hardcoded variables
	this_name = 'genfsys'
	indent = " "*len(header(this_name))
	version = 3.1
	other_versions_allowed = True

	# Run command line parser
	options, filters = main_parser()

	# Define folder names
	install_path = options['install_dir'] + '/'
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

	options_filename = locations['FSYS']['mainpath'] + '.options.dat'
	options_file = open(options_filename, 'w')
	for x in list(options.keys()):
		options_file.write("{0}\t\t{1}\n".format(x, options[x]))
	options_file.close()

	return options, filters, locations


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
			fields = line.split()
			if fields[0] not in locations:
				locations[fields[0]] = {}
			locations[fields[0]][fields[1]] = fields[2]
	return locations
