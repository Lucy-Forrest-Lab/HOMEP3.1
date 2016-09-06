# Name: genclib.py
# Language: python3
# Libraries:
# Description: Generates HOMEP chain library
# Author: Edoardo Sarti
# Date: Aug 15 2016

import os, shutil
from support import *

def from3to1(resname):
	f3t1 = {'ALA' : 'A',
	        'ARG' : 'R',
	        'ASN' : 'N',
	        'ASP' : 'D',
	        'CYS' : 'C',
	        'GLN' : 'Q',
	        'GLU' : 'E',
	        'GLY' : 'G',
	        'HIS' : 'H',
	        'ILE' : 'I',
	        'LEU' : 'L',
	        'LYS' : 'K',
	        'MET' : 'M',
	        'PHE' : 'F',
	        'PRO' : 'P',
	        'SER' : 'S',
	        'THR' : 'T',
	        'TRP' : 'W',
	        'TYR' : 'Y',
	        'VAL' : 'V'}

	if resname in list(f3t1.keys()):
		return f3t1[resname]
	else:
		return '0'

# PDB parser
def PDB_parser(locations, struct):
	struct_filename = locations['FSYS']['mainpath'] + locations['FSYS']['rpdb'] + struct + '.pdb'
	struct_file = open(struct_filename, 'r')
	text = struct_file.read().split("\n")
	struct_file.close()

	PDB_dict = {}
	res_ids = {}
	b_factor = {}
	b_norm = {}
	tech_list = ['NMR', 'X-RAY', 'THEORETICAL', 'ELECTRON']
	for line in text:
		if line[0:6] == 'EXPDTA':
			for tech in tech_list:
				if tech in line:
					PDB_dict['TECHNIQUE'] = tech
			if 'TECHNIQUE' not in PDB_dict:
				PDB_dict['TECHNIQUE'] = 'OTHER'
		elif line[0:10] == 'REMARK   2' and 'RESOLUTION' in line:
			fields = line.split()
			for nf in range(len(fields)):
				if 'ANGSTROM' in fields[nf]:
					if fields[nf-1] != "NULL":
						PDB_dict['RESOLUTION'] = float(fields[nf-1])
					else:
						PDB_dict['RESOLUTION'] = 0
			if 'RESOLUTION' not in PDB_dict and PDB_dict['TECHNIQUE'] != 'THEORETICAL' and PDB_dict['TECHNIQUE'] != 'NMR':
				raise NameError("ERROR: Resolution annotation of pdb {0} is badly formatted: {1}".format(struct_filename, line))
		elif line[0:10] == 'REMARK   3' and 'FREE R VALUE' in line and 'ERROR' not in line and 'SET' not in line and (line.split()[3] == 'FREE' or line.split()[3] == 'BIN'):
			try:
				PDB_dict['RFACTOR'] = float(line.split()[-1])
			except ValueError:
				PDB_dict['RFACTOR'] = 'NULL'
		elif line[0:4] == 'ATOM':
			if not line[21]:
				raise NameError("ERROR: There is an ATOM without chain name: {0}".format(line))
			ch_name = line[21].upper()
			if ch_name not in res_ids:
				res_ids[ch_name] = []
				b_factor[ch_name] = 0
				b_norm[ch_name] = 0
			if ch_name in res_ids and (not res_ids[ch_name] or  res_ids[ch_name][-1][0] != int(line[22:26])):
				res_ids[ch_name].append((int(line[22:26]), line[17:20]))
			if line[60:66]:
				b_factor[ch_name] += float(line[60:66])
				b_norm[ch_name] += 1
		elif line[0:5] == 'TITLE':
			if 'TITLE' not in PDB_dict:
				PDB_dict['TITLE'] = line[10:].rstrip()
			else:
				PDB_dict['TITLE'] += line[10:].rstrip()

	if 'RFACTOR' not in PDB_dict:
		PDB_dict['RFACTOR'] == 'NULL'

	PDB_dict['CHAIN'] = {}
	for chain in list(res_ids.keys()):
		if b_factor[chain] == 0:
			b_factor[chain] = 'NULL'
		else:
			b_factor[chain] = b_factor[chain] / b_norm[chain]
		PDB_dict['CHAIN'][chain] = [{'CHAINID' : chain,
	                                     'AVG_BFACTOR' : b_factor[chain],
		                             'NRES' : len(res_ids[chain])},
		                            {'RESIDS' : [x[0] for x in res_ids[chain]],
		                             'RESNAMES' : [from3to1(x[1]) for x in res_ids[chain]]}]

	return PDB_dict


# Structure checker
def checker(locations, database, filters):
	instructions = {}
	instructions_filename = locations['FSYS']['mainpath'] + '.superfamily_classification.dat'
	instructions_file = open(instructions_filename, 'w')
	exclusions_filename = locations['FSYS']['mainpath'] + locations['FSYS']['cpdb'] + 'exclusions.txt'
	exclusions_file = open(exclusions_filename, 'w')
	tab_filename = locations['FSYS']['mainpath'] + locations['FSYS']['cpdb'] + 'info.txt'
	tab_file = open(tab_filename, 'w')
	tab_string = ""
	new_database = {}
	for struct in list(database.keys()):
		PDB_dict = PDB_parser(locations, struct)

		for chain in PDB_dict['CHAIN']:
			exclude_chain = []
#			print(struct, chain)
#			print(list(database[struct][1]['CHAIN'].keys()))

			# Not present in PDB_TM oligomeric complex
			if chain not in database[struct][1]['CHAIN']:
				exclude_chain.append('Not contained in PDB_TM biological form')
			else:
				n_pdbtm = int(database[struct][1]['CHAIN'][chain][0]['NUM_TM'])
				# Is it alpha or beta?
#				print(database[struct][1]['CHAIN'][chain][0]['TYPE'])
				s_type = database[struct][1]['CHAIN'][chain][0]['TYPE']
				if s_type != 'alpha' and s_type != 'beta':
					exclude_chain.append('Chain does not cross the membrane')
			
			# NMR check
			if not filters['NMR'] and PDB_dict['TECHNIQUE'] == 'NMR':
				exclude_chain.append('NMR structure')

			# Theoretical model check
			if not filters['THM'] and PDB_dict['TECHNIQUE'] == 'THEORETICAL':
				exclude_chain.append('Theoretical model')

			# Holes greater than hole threshold check
			if filters['hole_thr'] > 0:
				for n in range(1,len(PDB_dict['CHAIN'][chain][1]['RESIDS'])):
					if PDB_dict['CHAIN'][chain][1]['RESIDS'][n] - PDB_dict['CHAIN'][chain][1]['RESIDS'][n-1] > filters['hole_thr']:
						exclude_chain.append('Contains hole longer than {0} residues'.format(filters['hole_thr']))
						break

			# Strange residues
			for res in PDB_dict['CHAIN'][chain][1]['RESNAMES']:
				if res == '0':
					exclude_chain.append('One or more residues with unsupported names')
					break

			# Resolution check
			if 'resolution' in filters and PDB_dict['TECHNIQUE'] != 'THEORETICAL' and PDB_dict['TECHNIQUE'] != 'NMR':
				if PDB_dict['RESOLUTION'] > filters['resolution']:
					exclude_chain.append('Resolution is higher than {0}'.format(filters['resolution']))
				elif PDB_dict['RESOLUTION'] == 0:
					exclude_chain.append('No resolution information found'.format(filters['resolution']))

			# Was there something wrong?
			if not exclude_chain:
				s_type = database[struct][1]['CHAIN'][chain][0]['TYPE']
				if struct not in instructions:
					instructions[struct] = {}
				instructions[struct][chain] = (s_type, n_pdbtm)
				instructions_file.write("{0}\t{1}\t{2}\n".format(s_type, n_pdbtm, struct+'_'+chain))
				chain_filename = locations['FSYS']['mainpath'] + locations['FSYS']['cpdb'] + struct + '_' + chain + '.pdb'
				chain_file = open(chain_filename, 'w')
				struct_filename = locations['FSYS']['mainpath'] + locations['FSYS']['rpdb'] + struct + '.pdb'
				struct_file = open(struct_filename, 'r')
				text = struct_file.read().split('\n')
				struct_file.close()
				for line in text:
					if line[0:4] == 'ATOM' and line[21] == chain:
						chain_file.write(line + '\n')
				chain_file.close()
				if struct not in new_database:
					new_database[struct] = []
					new_database[struct].append(database[struct][0])
					new_database[struct].append({})
					for key in [x for x in list(database[struct][1].keys()) if x != 'CHAIN']:
						new_database[struct][1][key] = database[struct][1][key]
					new_database[struct][1]['CHAIN'] = {}
				new_database[struct][1]['CHAIN'][chain] = database[struct][1]['CHAIN'][chain]
					
			else:
#				print(exclude_chain)
				exclusions_file.write(struct + '_' + chain + '\t\t' + exclude_chain[0] + '\n')
				for nl in range(1, len(exclude_chain)):
					exclusions_file.write(' '*len(struct) + ' ' + ' '*len(chain) + '\t\t' + exclude_chain[nl] + '\n')

		# Introduce PDB_dict as a key of the database
		if struct in new_database:
			new_database[struct][1]['FROM_PDB'] = PDB_dict

			if not tab_string:
				for key in [x for x in sorted(list(new_database[struct][1]['FROM_PDB'].keys())) if x != 'CHAIN']:
					if key == 'TITLE':
						tab_string += "{0:150} ".format(key)
					else:
						tab_string += "{0:16} ".format(key)
				tab_string += "{0:5} ".format('CHAIN')
				for key in [x for x in sorted(list(new_database[struct][1]['FROM_PDB']['CHAIN'][chain][0].keys())) if x != 'CHAINID']:
					tab_string += "{0:16} ".format(key)
				tab_file.write(tab_string[:-1] + "\n")

			tab_string = ""
			for key in [x for  x in sorted(list(new_database[struct][1]['FROM_PDB'].keys())) if x != 'CHAIN']:
				tab_string += "{0:16} ".format(new_database[struct][1]['FROM_PDB'][key])
				if key == 'TITLE':
					tab_string += "{0:150} ".format(new_database[struct][1]['FROM_PDB'][key][:150])
				else:
					tab_string += "{0:16} ".format(str(new_database[struct][1]['FROM_PDB'][key]))
			counter = 0
			strlen = len(tab_string)
			for chain in sorted(list(new_database[struct][1]['FROM_PDB']['CHAIN'])):
				if counter > 0:
					tab_string = " "*strlen
				tab_string += "{0:5} ".format(str(new_database[struct][1]['FROM_PDB']['CHAIN'][chain][0]['CHAINID']))
				for key in [x for x in sorted(list(new_database[struct][1]['FROM_PDB']['CHAIN'][chain][0].keys())) if x != 'CHAINID']:
					tab_string += "{0:16} ".format(str(new_database[struct][1]['FROM_PDB']['CHAIN'][chain][0][key]))
				tab_file.write(tab_string[:-1] + "\n")
				counter += 1
			tab_string = "XXX"

	exclusions_file.close()
	instructions_file.close()
	tab_file.close()

	return new_database, instructions


def structure_sorter(locations, instructions):
	ssd = {'alpha' : 'a', 'beta' : 'b'}
	for struct in instructions:
		for chain in instructions[struct]:
			ss = instructions[struct][chain][0]
			ntm = instructions[struct][chain][1]
			destination_dir = locations['FSYS']['mainpath'] + ss + '/' + str(ntm) + '/'
			if not os.path.exists(destination_dir):
				os.mkdir(destination_dir)
				os.mkdir(destination_dir + 'structures/')
			shutil.copy(locations['FSYS']['mainpath'] + locations['FSYS']['cpdb'] + struct + '_' + chain + '.pdb', destination_dir + 'structures/')


# Library function
def generate_chain_pdb_files(filters, locations, database):
	# Hardcoded variables
	this_name = 'genclib'
	indent = " "*len(header(this_name))
	version = 3.1

	# Checks
	for path_name in [locations['FSYS']['mainpath'] + locations['FSYS'][x] for x in list(locations['FSYS'].keys()) if x != 'installpath' and x != 'mainpath' and x!= 'main']:
		if not os.path.exists(path_name):
			logmsg = header(this_name) + "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name)
			write_log(this_name, logmsg)	
			raise NameError(logmsg)

	# Structure checker
	database, fs_ordering = checker(locations, database, filters)

	# Filesystem sorting
	structure_sorter(locations, fs_ordering)

	return database
