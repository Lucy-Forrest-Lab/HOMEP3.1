# Name: generate_library.py
# Language: python3
# Libraries: argparse, genfsys, genrlib, genclib, straln, clusterize
# Description: Updates HOMEP library
# Author: Edoardo Sarti
# Date: Aug 17 2016

#import os, sys, multiprocessing, subprocess, re, time

import argparse, genfsys, genrlib, genclib, straln, clusterize

# parser and checks

# python generate_library.py -d -pdbtm -rf 3.5
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--main_dir', nargs=1)
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

#Add check if main_dir path exists


filters = {'resolution' : float(parsed.resolution_filter[0]),
           'NMR' : parsed.with_nmr,
           'THM' : parsed.with_theoretical,
           'hole_thr' : int(parsed.hole_thr[0])}

# execute

locations = genfsys.filesystem_info(str(parsed.main_dir[0]))

pdbtm_data, diff_database_namelist = genrlib.update_raw_pdb_library(locations, str(parsed.pdbtm_file_path[0]))
# PDB names must be in upper case
# After downloading, check for existence and then compile a list and a no-list
# In pdbtm_data output, there must be all info regarding the structures

if not pdbtm_data:
	raise NameError("ERROR: pdbtm_data not produced")

if diff_database_namelist:
	diff_database = {}
	for struct in diff_database_namelist:
		diff_database[struct] = pdbtm_data[struct]

	diff_database = genclib.generate_chain_pdb_files(locations, diff_database, filters)
	# Here, operate any possible checks. The resulting list must be the cleanest possible
	# After checking, filter by resolution, then divide by number of TM domains, then create filesystem and add codes
	# Eventually there must be two folders: one with all identified chains, another with the used chains

	for struct in diff_database_namelist:
		pdbtm_data[struct] = diff_database[struct]

	np = int(parsed.number_of_procs[0])

	table = straln.structure_alignment(locations, str(parsed.straln_path[0]), np, str(parsed.output_tab[0]))
	# Take this part from start_FrTM.py and adapt
else:
	table = straln.make_new_table(locations, str(parsed.output_tab[0]))

homep_library = clusterize.clusterize(locations, pdbtm_data, table, str(parsed.output_tab[0]), str(parsed.HOMEP_filename[0]), float(parsed.object_thr[0]), float(parsed.cluster_thr[0]))
# Must report each used chain
# Must output library table
