# Name: generate_library.py
# Language: python3
# Libraries: genfsys, genrlib, genclib, straln, clusterize
# Description: Generates HOMEP library from scratch
# Author: Edoardo Sarti
# Date: Sep 4 2016

import genfsys
import genrlib
import genclib
import straln
import clusterize


print("GENERATE FILESYSTEM")
# Parse the command line and generate filesystem.
# The main directory of the filesystem is created in the installation
# path and is named 'HOMEP_3.1_YYYY_MM_DD', where YYYY, MM and DD are
# the year, month and day of creation.
# Return three lists:
#  options - contains all command line options
#  filters - contains all activated filters (see genclib)
#  locations - contains all addresses of the filesystem locations
# Write two hidden files in the main folder of the fileystem:
#  .options - transcription of the options list
#  .locations - transcription of the locations list
options, filters, locations = genfsys.generate_filesystem()


print("GENERATE RAW LIBRARY")
# Parse the PDTBM library file, download pdb files from PDB, builds a
# comprehensive nested structure containing all information from the PDBTM
# library file.
# PDB codes are stored in upper case (PDB-style).
# Return a nested structure:
#  pdbtm_data - contains all information parsed from the PDBTM library file
#               It is a dictionary whose keys are the PDB 4-letter codes
#               (capitalized). Each entry is a list of two elements, the first
#               containing the dictionary of general information keywords (the
#               ones contained in the opening tag) and the second containing
#               the dictionary of specific information tags (the nested tags).
#               The latter is build recursively in the same way (it is itself
#               a two-element list). When a tag is non-unique (i.e., there can
#               be more than one 'CHAIN' tag), the dictionary entry contains a
#               list, with one element for each tag. Each element is again a
#               two-element list...
pdbtm_data = genrlib.generate_raw_pdb_library(options, locations)


print("GENERATE CHAIN LIBRARY")
# Thoroghly check the pdb files, divide them into chains, selects chains
# according to the filters, compiles a smaller database containing only the
# selected chains.
# The selected chains are copied into the pdb chain folder and into the
# proper superfamily folder (in its structure/ subfolder).
# Return the nested structure:
#  core_pdbtm_data - structured as pdbtm_data, but only containing the
#                    selected chains. In addition to the information from the
#                    PDBTM library file, each selected chain will have a new
#                    'FROM_PDB' keyword holding a dictionary with the
#                    information retrieved from the pdb file.
core_pdbtm_data = genclib.generate_chain_pdb_files(filters, locations, pdbtm_data)


print("CALCULATE STRUCTURE ALIGNMENTS")
table = straln.structure_alignment(options, locations)
# Take this part from start_FrTM.py and adapt


print("CLUSTERIZE RESULTS")
homep_library = clusterize.clusterize(options, locations, core_pdbtm_data, table)
# Must report each used chain
# Must output library table
