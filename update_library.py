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
#  .options.dat - transcription of the options list
#  .locations.dat - transcription of the locations list
options, filters, locations = genfsys.filesystem_info()


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
pdbtm_data, diff_pdbtm_data = genrlib.upate_raw_pdb_library(options,
                                                            locations)

exit(1)

print("GENERATE CHAIN LIBRARY")
# Thoroughly check the pdb files, divide them into chains, selects chains
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
# Write one hidden file in the main folder of the filesystem:
#  .superfamily_classification.dat - contains the proper location for each pdb
#                                    structure
# Write two log files in the chain pdb folder:
#  exclusions.txt - states the reasons why each excluded chain was ruled out
#  info.txt - contains an unclassified summary of all selected chains,
#             providing resolution, technique, rfactor, title and chain-wise
#             specifics for each element
core_pdbtm_data = genclib.generate_chain_pdb_files(filters,
                                                   locations,
                                                   pdbtm_data)


print("CALCULATE STRUCTURE ALIGNMENTS")
# Check the existent alignments, write the alignments to perform, runs the
# specified alignment program in parallel over the specified number of
# processors. Each job performs all due alignments with a specified structure
# For example, if the alignments A-B, B-A, C-A, C-B have to be run, the 
# jobs will be (A-B), (B-A), (C-A, C-B). When an alignment is complete, add
# to the corresponding sequence andn structure alignment files the new data.
# Files are organized in the alignments/ folder contained in each superfamily
# directory. In the alignments/ folder there is the fasta/ subfolder that 
# contains the sequence alignments, and the str_alns/ subfolder that contains
# the structural alignments. The alignments are organized in files carrying
# the name of the first structure in the alignment. For example, the sequence
# alignments A-B and A-H will be contained in the file seq_A.dat, while the
# structure alignment C-A will be contained in the file str_C.dat.
# Return the nested structure:
#  table - contains the sequence identity, TM-score and RMSD of each alignment
#          present in the database (not just the ones recently performed).
#          Its access method is: table[clade][superfamily][struct_1][struct_2]
# Write the file in the main folder of the file system:
#  <structure_alignments.dat> - the name of this file can be changed with the
#                               flag --output_tab. Contains the information
#                               placed in the nested structure 'table'
table = straln.structure_alignment(options,
                                   locations)
# Take this part from start_FrTM.py and adapt


print("CLUSTERIZE RESULTS")
# Starting from the data contained in the nested structure 'table' (or the
# corresponding file) divide the structures in each superfamily into families
# and groups.
# Return the nested structure:
#  homep_library - contains the complete HOMEP3.1 library containing all data
#                  present in the filesystem. The structures are organized in
#                  Clades (alpha, beta), Superfamilies (number of TM domains),
#                  Families (proximity with TM-score metric), Objects (proximity
#                  with sequence identity matrix).
# Write the file in the main folder of the file system:
#  <HOMEP3.1.dat> - the name of this file can be changed with the flag
#                   --output_homep. Contains the full classification of all
#                   structures contained in the database.
homep_library = clusterize.clusterize(options, locations, core_pdbtm_data, table)
