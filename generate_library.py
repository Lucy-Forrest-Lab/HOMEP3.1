# Name: generate_library.py
# Language: python3
# Libraries: argparse, genfsys, genrlib, genclib, straln, clusterize
# Description: Generates HOMEP library from scratch
# Author: Edoardo Sarti
# Date: Aug 10 2016

#import os, sys, multiprocessing, subprocess, re, time

import argparse, genfsys, genrlib, genclib, straln, clusterize

# parser and checks

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--install_dir', nargs=1)
parser.add_argument('-pdbtm', '--pdbtm_file_path', nargs=1)
parser.add_argument('-s', '--straln_path', nargs=1)
parser.add_argument('-np', '--number_of_procs', nargs=1)
parser.add_argument('-ot', '--object_thr', nargs=1)
parser.add_argument('-ct', '--cluster_thr', nargs=1)
parser.add_argument('-o', '--output_file', nargs=1)
np = int(number_of_procs)

#Add check if main_dir path exists

# execute

locations = genfsys.generate_filesystem(install_dir)

pdbtm_data = genrlib.generate_raw_pdb_library(locations, pdbtm_file_path)
# PDB names must be in upper case
# After downloading, check for existence and then compile a list and a no-list
# In pdbtm_data output, there must be all info regarding the structures

genclib.generate_chain_pdb_files(locations, pdbtm_data)
# Here, operate any possible check. The resulting list must be the cleanest possible
# After checking, filter by resolution, then divide by number of TM domains, then create filesystem and add codes
# Eventually there must be two folders: one with all identified chains, another with the used chains

straln.structure_alignment(locations, straln_path, np)
# Take this part from start_FrTM.py and adapt

homep_library = clusterize.clusterize(locations, pdbtm_data, object_thr, cluster_thr)
# Must report each used chain
# Must output library table
