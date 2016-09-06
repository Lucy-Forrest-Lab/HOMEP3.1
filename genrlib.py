# Name: genrlib.py
# Language: python3
# Libraries: re, urllib.request, gzip, shutil, os, copy, support 
# Description: Generates HOMEP raw pdb library
# Author: Edoardo Sarti
# Date: Sep 05 2016

import re
import urllib.request
import gzip
import shutil
import os
import copy
from support import *


# --- Support functions for XML parser ------------------------------------- #

# Parse attributes in an opening or opening/closing tag
def parse_attributes(line):
	attr = {}
	line = line.strip().replace('/>','').replace('>','')
	for field in line.split():
		if '=' in field:
			attr[field.split('=')[0]] = field.split('=')[1][1:-1]
	if attr:
		return attr
	else:
		return None


# Parse tag name in any type of tag
def extract_tag(line, this_name):
	# Extracts the first word in the first tag (opening or closing) found in the line
	tag = re.findall(r'<(.*?)/*>',line,re.DOTALL)
	if not tag:
		logmsg = header(this_name) + "ERROR: no tag name found in line {0}".format(line)
		write_log(this_name, logmsg)	
		raise NameError(logmsg)
	return tag[0].split()[0]


# Parse text outside tags
def extract_text(line, this_name):
	# Text between opening and closing tags
	text = re.findall(r'>(.*\S.*)</',line,re.DOTALL)
	if not text:
		# Text after an opening tag
		text = re.findall(r'>(.*\S.*)$',line,re.DOTALL)
	if not text:
		# Text before a closing tag
		text = re.findall(r'^(.*\S.*)</',line,re.DOTALL)
	if not text:
		logmsg = header(this_name) + "ERROR: no text found in line {0}".format(line)
		write_log(this_name, logmsg)
		raise NameError(logmsg)
	return text[0]

# --- End support functions for XML parser ---------------------------------- #


# XML parser
# Parse the PDBTM library file.
# Fill a nested dictionary structure composed in the following way:
#  - The main keys of the dictionary are the (capitalized) 4-letter PDB codes
#  - A tag is rendered as a list of two dictionaries. The former contains all keywords
#    that appear in the opening tag, the latter contains everything is nested into the
#    tag (including other tags, their keywords, their contents, and so on).
#  - When the tag is allowed to repeat more than once it must appear in the open_list_tags
#    hardcoded list. If a tag of this sort is encountered, the keyword of the dictionary
#    corresponds to a list of two-element lists - not directly to the two-element list of 
#    the tag.
#  - The only two exceptions to this set of rules are the opening tag (i.e., the main
#    keyword is not 'pdbtm', as it should, but are the 4-letter PDB codes) and the
#    'CHAIN' tag, which does not hold a list of chains as it should, rather contains a
#    dictionary whose keywords are the names of the chains. This makes future searches in
#    the database easier.
# Example: in order to know all the reported names of the chains of the structure 
# corresponding to the PDB code 1AAA, one should query database['1AAA'][1]['CHAIN'].keys()
# (Note that this structure is acquired only at the end of the process!)
def parser(pdbtm_file_path, this_name):
	# The resulting object is a dictionary where each 4-letter PDB name corresponds
	# to all information contained in the database. Each PDB name corresponds to a
	# 2-entry list, the first entry being the header, the second the body.
	# Data are nested following the same order of the database. For example, in order
	# to verify how many TM domains has the 3rd chain of the structure 1A0S, we have
	# to read the entry database['1A0S'][1]['CHAIN'][2]['NUM_TM'], because we want to
	# access the 'NUM_TM' key of the 3rd item of the list corresponding to the 'CHAIN'
	# key of the body (1) dictionary of structure 1A0S.
	# The database is compiled in a slightly cumbersome yet recursive manner, which is
	# explained case by case in the code. Then, it is translated into the final structure.
	pdbtm_file = open(pdbtm_file_path, 'r')
	text = pdbtm_file.read().split('\n')
	pdbtm_file.close()

	DB = [[None, {}]]
	DB_tagnames = []
	open_list_tags = ['pdbtm', 'MODIFICATION', 'MATRIX', 'CHAIN', 'REGION']
	stophere = 0

	for line in text:
		if not line:
			continue
#		if '<pdbtm' in line and DB and 'pdbtm' in DB[1][1]:
#			print(DB[1][1]['pdbtm'][-1], DB_tagnames)
		# If line contains a tag
		if re.search('<.*>', line):
			# If line is a comment, skip
			if re.search('^\s*<\?', line):
#				print("comment:\t"+line)
				continue
			# If line contains a tag but does not contain a closing tag (i.e. only contains one or more opening tags)
			elif not re.search('</.*>', line):
				# If line is a standalone tag, compile a "body termination" dictionary
				# containing the parameters inside the standalone tag, and extract
				# the name of the tag.
				if re.search('^\s*<.*/>\s*$', line):
#					print("O/C:\t"+line)
					parameters = parse_attributes(line)
					tag = extract_tag(line, this_name)
					# If the tag accepts multiple instances
					if tag in open_list_tags:
						# If the tag is not present yet in the body dictionary
						# of the last element of the DB list, place a new element
						# whose key is the tag name and whose value is an empty list.
						if tag not in DB[-1][1]:
							DB[-1][1][tag] = []
						# In the body dictionary of the last element of the DB list,
						# append to the list corresponding to the tag name the "body
						# termination" dictionary
						DB[-1][1][tag].append(parameters)
					# If the tag does not accept multiple instances, just add an element to
					# the body dictionary of the last element of the DB list. The key of this
					# element is the tag name and the value is the "body termination" dictionary
					else:
						DB[-1][1][tag] = parameters
				# If line is an unmatched opening tag, compile a header dictionary 
				# containing the parameters inside the opening tag, and append to the
				# DB list a new element composed of the header and the empty body,
				# shaped as a dictionary (it will be changed in case the body is text).
				# In the DB_tagnames list, take note of the name of the unclosed tag.
				elif re.search('^\s*<.*>\s*$', line):
#					print("O:\t"+line)
					parameters = parse_attributes(line)
					element = [parameters, {}]
					tag = extract_tag(line, this_name)
					DB.append(element)
					DB_tagnames.append(tag)
				# If line contains an opening tag at the beginning
				elif re.search('^\s*<.*>.*$', line):
#					print("OT:\t"+line)
					tag = extract_tag(line, this_name)
					text = extract_text(line, this_name).strip() + " "
					element = [None, text]
					DB.append(element)
					DB_tagnames.append(tag)
				else:
					raise NameError("ERROR: line {0} is not compliant with the XML format".format(line))
			# If line contains a tag and contains a closing tag
			else:
				# If line contains opening and closing tags, extract the tag and the text,
				# and, in the body dictionary of the last element of the DB list, add the 
				# text as the value corresponding to the key tag.
				if re.search('^\s*<.*>.+</.*>\s*$', line):
#					print("O+C:\t"+line)
					tag = extract_tag(line, this_name)
					text = extract_text(line, this_name)
					DB[-1][1][tag] = text
				# If line does not contain opening and closing tags
				else:
					# If line is a closing tag
					if re.search('^\s*</.*>\s*$', line):
#						print("C:\t"+line)
						pass
					# If line contains a closing tag at the end
					elif re.search('^.+</.*>\s*$', line):
#						print("TC:\t"+line)
						# If there is no accumulated text yet and the body of the last element of the DB list is still an empty
						# dictionary, change it to text and add the first line.
						if type(DB[-1][1]) == dict:
							DB[-1][1] = extract_text(line, this_name).strip() + " "
						# If there already is accumulated text, just extend the string
						else:
							DB[-1][1] += extract_text(line, this_name).strip() + " "
					else:
						raise NameError("ERROR: line {0} is not compliant with the XML format".format(line))
					# If the last tag name in the DB_tagnames list is a tag accepting multiple instances
					if DB_tagnames[-1] in open_list_tags:
						# If the tag is not present yet in the body dictionary of the second-last element
						# of the DB list, place a new element whose key is the tag name and whose value is 
						# an empty list
						if DB_tagnames[-1] not in DB[-2][1]:
							DB[-2][1][DB_tagnames[-1]] = []
						# In the body dictionary of the second-last element of the DB list, append to the
						# list corresponding to the tag name the whole last element of the DB list
						DB[-2][1][DB_tagnames[-1]].append(DB[-1])
					# If the last tag name in the DB_tagnames list is not a tag accepting multiple instances,
					# just add an element to the body dictionary of the second-last element of the DB list.
					# The key of this element is the tag name and the value is the whole last element of the DB list.
					else:
						DB[-2][1][DB_tagnames[-1]] = DB[-1]
					# Delete the last element of the DB and DB_tagnames lists, since now is incorporated as an element
					# of the body dictionary of the second-last element of the DB list.
					del DB_tagnames[-1], DB[-1]

		# If line does not contain tags
		else:
#			print("T:\t"+line)
			# If there is no accumulated text yet and the body of the last element of the DB list is still an empty
			# dictionary, change it to text and add the first line.
			if type(DB[-1][1]) == dict:
				DB[-1][1] = line.strip() + " "
			# If there already is accumulated text, just extend the string
			else:
				DB[-1][1] += line.strip() + " "

		# If you want to try with a finite number of structures, this command will only parse the first N structures of the
		# database file.
		'''
		N = 2
		if '</pdbtm>' in line:
			print(DB)
			stophere += 1
		if stophere == N:
			break
		'''
#		if '</PDBTM>' in line:
#			print("READING COMPLETE!")

	
	# Reorganize the database in a dictionary structure
	# The DB structure is a list with only one element (this is the result of recursion), whose body contains a dictionary with only
	# one element ('pdbtm') with N instances, where N is the number of structures in the database. Here, data is reorganized in a
	# dictionary having as keys the uppercase 4-letter PDB names, and for values the corresponding entry compiled in the DB[0][1]['pdbtm']
	# sublist.
	DB2 = {}
	for ns in range(len(DB[0][1]['PDBTM'][1]['pdbtm'])):
		pdbname = DB[0][1]['PDBTM'][1]['pdbtm'][ns][0]['ID'].upper()
		tmp_struct = copy.deepcopy(DB[0][1]['PDBTM'][1]['pdbtm'][ns])
		tmp_struct[1]['CHAIN'] = {}
		for nc in range(len(DB[0][1]['PDBTM'][1]['pdbtm'][ns][1]['CHAIN'])):
			chain = DB[0][1]['PDBTM'][1]['pdbtm'][ns][1]['CHAIN'][nc][0]['CHAINID'].upper()
			tmp_chain =  DB[0][1]['PDBTM'][1]['pdbtm'][ns][1]['CHAIN'][nc]
			tmp_struct[1]['CHAIN'][chain] = tmp_chain
		DB2[pdbname] = tmp_struct
	return DB2


def mini_parser(pdbtm_file_path, this_name):
	pdbtm_file = open(pdbtm_file_path, 'r')
	text = pdbtm_file.read().split('\n')
	pdbtm_file.close()

	DB = {}
	for line in text:
		if not line:
			continue
		if re.search('<pdbtm.*>', line):
			attr = parse_attributes(line)
			DB[attr['ID']] = [attr, {}]
	return DB


# Download structures from the PDB website
def download_structures(database_namelist, raw_pdb_dir):
	for pdbname in database_namelist:
		url = 'http://www.rcsb.org/pdb/files/'+pdbname+'.pdb.gz'
		local_filename = raw_pdb_dir + pdbname + '.pdb'
		if not os.path.exists(local_filename):
			with urllib.request.urlopen(url) as response:
				with gzip.GzipFile(fileobj=response) as uncompressed, open(local_filename, 'wb') as local_file:
					shutil.copyfileobj(uncompressed, local_file)

	downloaded_files = [x[:-4] for x in os.listdir(raw_pdb_dir) if x[-4:]=='.pdb']
	missing_files_filename = raw_pdb_dir + 'missing_files.txt'
	missing_files_file = open(missing_files_filename, 'a')
	for database_struct in database_namelist:
		if database_struct not in downloaded_files:
			missing_files_file.write(database_struct+"\n")
	missing_files_file.close()


# --- Library functions ---------------------------------------------------- #

def generate_raw_pdb_library(options, locations):
	# Hardcoded variables
	this_name = 'generate_raw_pdb_library'
	version = 3.1

	pdbtm_file_path = options['pdbtm_file_path']

	# Checks
	for path_name in [locations['FSYS']['mainpath']+locations['FSYS'][x] for x in list(locations['FSYS'].keys()) if x != 'installpath' and x != 'mainpath' and x!= 'main']:
		if not os.path.exists(path_name):
			raise_error(this_name, "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name))

	# Parser
	database = parser(pdbtm_file_path, this_name)
	database_namelist = list(database.keys())
	shutil.copy(pdbtm_file_path, locations['FSYS']['mainpath']+'pdbtm_database.dat')
	

	# Downloader
	download_structures(database_namelist, locations['FSYS']['mainpath'] + locations['FSYS']['rpdb'])

	return database

def update_raw_pdb_library(locations, pdbtm_file_path):
	# Hardcoded variables
	this_name = 'update_raw_pdb_library'
	version = 3.1
	log = ""

	# Checks
	for path_name in [locations['FSYS']['mainpath']+locations['FSYS'][x] for x in list(locations['FSYS'].keys()) if x != 'installpath' and x != 'mainpath' and x!= 'main']:
		if not os.path.exists(path_name):
			raise_error(this_name, "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name))

	# Parser
	# Parse new database
	database = parser(pdbtm_file_path, this_name)
	database_namelist = list(database.keys())

	# Parse old database and calculate the difference
	old_pdbtm_file_path = locations['FSYS']['mainpath']+'pdbtm_database.dat'
	old_database = parser(old_pdbtm_file_path, this_name)
	old_database_namelist = list(old_database.keys())

	# Overwrite the pdbtm_databse.dat file and save the old one as a hidden file
	# (until another update won't overwrite on it as well)
	shutil.copy(locations['FSYS']['mainpath']+'pdbtm_database.dat', locations['FSYS']['mainpath']+'.previous_pdbtm_database.dat')
	shutil.copy(pdbtm_file_path, locations['FSYS']['mainpath']+'pdbtm_database.dat')

	diff_database_namelist = []
	for struct in database_namelist:
		if struct not in old_database_namelist:
			diff_database_namelist.append(struct)

	if not diff_database_namelist:
		log += print_log(this_name, "No updates to be done.")
		return database, []
	
	# Downloader
	download_structures(diff_database_namelist, locations['FSYS']['mainpath'] + locations['FSYS']['rpdb'])

	return database, diff_database_namelist
