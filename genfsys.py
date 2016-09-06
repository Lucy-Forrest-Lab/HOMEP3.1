# Name: genfsys.py
# Version: 3.1
# Language: python3
# Libraries: sys, os, datetime, argparse, support
# Description: Generates HOMEP file system
# Author: Edoardo Sarti
# Date: Sep 05 2016

import sys
import os
import datetime
import argparse
from support import *


# Command line parser
# Return two dictionaries:
#  options - contains all options specified in the command line and/or in 
#            defaults
#  filters - contains all optional or adjustable criteria of selection that
#            will be used to filter the database in genclib.
def main_parser():
	this_name = 'main_parser'
	parser = argparse.ArgumentParser()

	# Compulsory arguments
	parser.add_argument('-pdbtm', '--pdbtm_file_path', nargs=1)

	# Either of these must be set
	parser.add_argument('-d', '--install_dir', nargs=1)  # Either of these
	parser.add_argument('-m', '--main_dir', nargs='?')   # must be set

	# These are compulsory if -d is set, will not be considered if -m is set
	parser.add_argument('-s', '--straln_path', nargs='?')
	parser.add_argument('-np', '--number_of_procs', nargs='?')
	parser.add_argument('-ot', '--object_thr', nargs='?')
	parser.add_argument('-ct', '--cluster_thr', nargs='?')
	parser.add_argument('-rf', '--resolution_filter', nargs='?')

	# Optional
	parser.add_argument('-ht', '--hole_thr', nargs='?')
	parser.add_argument('-oh', '--output_homep', nargs='?')
	parser.add_argument('-otab', '--output_tab', nargs='?')
	parser.add_argument('-with_nmr', action='store_true')
	parser.add_argument('-with_theoretical', action='store_true')

	# Default values for optional arguments
	parser.set_defaults(install_dir = '')
	parser.set_defaults(main_dir = '')
	parser.set_defaults(hole_thr = '100')
	parser.set_defaults(output_tab = 'structure_alignments.dat')
	parser.set_defaults(output_homep = 'HOMEP3.1.dat')

	parsed = parser.parse_args()


	# Create 'options' dictionary containing all selected options as
	# strings
	options = {}
	for x in sorted(parsed.__dict__):
		if type(parsed.__dict__[x]) == list:
			print(x, parsed.__dict__[x][0])
			options[x] = parsed.__dict__[x][0]
		else:
			print(x, parsed.__dict__[x])
			options[x] = parsed.__dict__[x]

	for x in list(options.keys()):
		if options[x] == 'False':
			options[x] = False
		elif options[x] == 'True':
			options[x] = True

	# Either flag -d or -m must be set
	if not (options['install_dir'] or options['main_dir']):
		raise_error(this_name, "ERROR: either -d (--install_dir) or -m (--main_dir) must be specified")

	# If -d flag is set, these flags must be set too
	if options['install_dir'] and not (options['straln_path'] or options['number_of_procs'] or 
	 options['object_thr'] or options['cluster_thr'] or options['resolution_filter']):
		raise_error(this_name, "ERROR: if flag -d (--install_dir) is set, these flags must be set too:\n" +
		                       "       -s  (--straln_path) <structure_aln_program_path>\n" +
		                       "       -np (--number_of_procs) <number of processors>\n" +
		                       "       -ot (--object_thr) <seq id threshold to determine Object>\n" +
		                       "       -ct (--cluster_thr) <TM-score threshold to clusterize Objects>\n" +
		                       "       -rf (--resolution_filter) <resolution filter (in A) to select chains>\n")

	# Create 'filters' dictionary containing all optional and tunable
	# criteria of selection for genclib
	filters = {}
	if 'resolution_filter' in options and options['resolution_filter']:
		filters['resolution'] = float(options['resolution_filter'])
	if 'with_nmr' in options and options['with_nmr']:
		filters['NMR'] = options['with_nmr']
	if 'with_theoretical' in options and options['with_theoretical']:
		filters['THM'] = options['with_theoretical']
	if 'hole_thr' in options and options['hole_thr']:
		filters['hole_thr'] = int(options['hole_thr'])

	return options, filters


def write_hidden_files(options, filters, locations):
	# Write the .options.dat hidden file in the main directory of the file system.
	# It contains all command line options with whom the script is running.
	options_filename = locations['FSYS']['mainpath'] + '.options.dat'
	options_file = open(options_filename, 'w')
	for x in list(options.keys()):
		options_file.write("{0}\t\t{1}\n".format(x, options[x]))
	options_file.close()

	# Write the .filters.dat hidden file in the main directory of the file system.
	# It contains all filters with whom the script is running.
	filters_filename = locations['FSYS']['mainpath'] + '.filters.dat'
	filters_file = open(filters_filename, 'w')
	for x in list(filters.keys()):
		filters_file.write("{0}\t\t{1}\n".format(x, filters[x]))
	filters_file.close()

	# Write the .locations.dat hidden file in the main directory of the file system.
	# It contains the same information contained in the 'locations' dictionary.
	locations_filename = locations['FSYS']['mainpath'] + '.locations.dat'
	locations_file = open(locations_filename, 'w')
	for x in list(locations['FSYS'].keys()):
		locations_file.write("{0}\t\t{1}\t\t{2}\n".format('FSYS', x, locations['FSYS'][x]))
	locations_file.close()

	return

def read_and_merge(filenames, dictionaries):
	options_filename, filters_filename, locations_filename = filenames
	options, filters, locations = dictionaries
	# Read the .options.dat file and returns the 'options' dictionary
	stored_options = {}
	options_file = open(options_filename, 'r')
	text = options_file.read().split('\n')
	for line in text:
		if line:
			fields = line.split()
			if len(fields) == 2:
				stored_options[fields[0]] = fields[1]

	for x in list(stored_options.keys()):
		if x not in options:
			options[x] = stored_options[x]

	for x in list(options.keys()):
        	if options[x] == 'False':
        		options[x] = False
        	elif options[x] == 'True':
        		options[x] = True

	# Read the .filters.dat file and returns the 'options' dictionary
	stored_filters = {}
	filters_file = open(filters_filename, 'r')
	text = filters_file.read().split('\n')
	for line in text:
		if line:
			fields = line.split()
			stored_filters[fields[0]] = fields[1]

	for x in list(stored_filters.keys()):
		if x not in filters:
			filters[x] = stored_filters[x]

	for x in list(filters.keys()):
		if filters[x] == 'False':
			filters[x] = False
		elif filters[x] == 'True':
			filters[x] = True

	# Read the .locations.dat file and returns the 'locations' dictionary
	locations = {}
	locations_file = open(locations_filename, 'r')
	text = locations_file.read().split('\n')
	for line in text:
		if line:
			fields = line.split()
			if fields[0] not in locations:
				locations[fields[0]] = {}
			locations[fields[0]][fields[1]] = fields[2]

	return options, filters, locations


# --- Library functions ---------------------------------------------------- #

# Generate the file system for the HOMEP database 
def generate_filesystem():
	# Hardcoded variables
	this_name = 'generate_filesystem'
	indent = " "*len(header(this_name))
	version = 3.1
	other_versions_allowed = True

	# Run command line parser
	options, filters = main_parser(this_name)

	# Define folder names
	install_path = options['install_dir'] + '/'
	main_dir = 'HOMEP_' + str(version) + '_' + datetime.datetime.now().strftime("%Y_%m_%d") + '/'
	main_path = install_path + main_dir
	rpdb_dir = 'raw_pdbs/'
	cpdb_dir = 'pdbs/'
	lib_dir = {}
	lib_dir['alpha'] = 'alpha/'
	lib_dir['beta'] = 'beta/'

	# Run checks over names and addresses
	if not os.path.exists(install_path):
		raise_error(this_name, "ERROR: The installation directory path {0} does not exist. Please specify an existing path.".format(install_path))
	if os.path.exists(main_path):
		raise_error(this_name, "ERROR: In the installation directory path {0} there already is a folder named {1}\n".format(install_path, main_dir))
	for filename in os.listdir(install_path):
		if filename[0:5] == 'HOMEP' and not other_versions_allowed:
			raise_error(this_name, "ERROR: In the installation directory path {0} there are other versions of HOMEP.\n".format(install_path) +
			              indent + "       If you want to continue, you have to set the internal variable other_versions_allowed as True.")

	# Generate filesystem
	log = ""
	os.mkdir(main_path)
	log += print_log(this_name, "Main directory created: {0}\n".format(main_path))

	os.mkdir(main_path + rpdb_dir)
	log += print_log(this_name, "Directory to store raw pdbs created: {0}\n".format(main_path + rpdb_dir))

	os.mkdir(main_path + cpdb_dir)
	log += print_log(this_name, "Directory to store curated pdbs created: {0}\n".format(main_path + cpdb_dir))

	os.mkdir(main_path + '.old/')
	log += print_log(this_name, "Hidden directory to store previous system files: {0}\n".format(main_path + '.old/'))

	for ss in 'alpha', 'beta':
		os.mkdir(main_path + lib_dir[ss])
		log += print_log(this_name, "Directory to store " + ss + " superfamilies created: {0}\n".format(main_path + lib_dir[ss]))

	write_log(this_name, log)

	# Create 'locations' nested dictionary.
	# Under the keyword 'FSYS' should go all paths and names relative to the file system;
	# Under the keyword 'OPT' should go all other locations and paths it's convenient to save.
	locations = {'FSYS' : {}, 'OPT' : {}}
	locations['FSYS']['installpath'] = install_path
	locations['FSYS']['mainpath'] = main_path
	locations['FSYS']['main'] = main_dir
	locations['FSYS']['rpdb'] = rpdb_dir
	locations['FSYS']['cpdb'] = cpdb_dir
	locations['FSYS']['alpha'] = lib_dir['alpha']
	locations['FSYS']['beta'] = lib_dir['beta']

	# Write the three hidden files	
	write_hidden_files(options, filters, locations)

	return options, filters, locations


# Retrieves all infromation about the file system

def filesystem_info():
	# Hardcoded variables
	this_name = 'filesystem_info'
	indent = " "*len(header(this_name))
	version = 3.1

	# Run command line parser
	options, filters = main_parser()
	main_path = options['main_dir']
	
	# Perform checks on paths
	if not os.path.exists(main_path):
		raise_error(this_name, "ERROR: Main directory {0} not found.".format(main_path))

	# If the .filters.dat, .locations.dat or the .options.dat files are not found, it returns error
	locations_filename = main_path + '/' + '.locations.dat'
	if not os.path.exists(locations_filename):
		raise_error(this_name, "ERROR: File {0} not found. Filesystem corrupted.".format(locations_filename))
	options_filename = main_path + '/' + '.options.dat'
	if not os.path.exists(options_filename):
		raise_error(this_name, "ERROR: File {0} not found. Filesystem corrupted.".format(options_filename))
	filters_filename = main_path + '/' + '.filters.dat'
	if not os.path.exists(filters_filename):
		raise_error(this_name, "ERROR: File {0} not found. Filesystem corrupted.".format(filters_filename))

	locations = {}
	options, filters, locations = read_and_merge([options_filename, filters_filename, locations_filename], [options, filters, locations])
	archive_old_file(locations, [options_filename, filters_filename, locations_filename])
	write_hidden_files(options, filters, locations)

	return options, filters, locations
