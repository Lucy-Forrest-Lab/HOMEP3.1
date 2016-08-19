# Name: genclib.py
# Language: python3
# Libraries:
# Description: Generates HOMEP chain library
# Author: Edoardo Sarti
# Date: Aug 15 2016

import ppm_segments, os, shutil 
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


# Structure checker
def checker(locations, database, filters):
	opm_data = {}
	if filters['OPM_TMdoms']:
		opm_data = ppm_segments.OPM_TMdoms()

#	print(opm_data)
#	raise NameError("HHH")

	instructions = {}
	exclusions_filename = locations['mainpath'] + locations['cpdb'] + 'exclusions.txt'
	exclusions_file = open(exclusions_filename, 'w')
	tab_filename = locations['mainpath'] + locations['cpdb'] + 'info.txt'
	tab_file = open(tab_filename, 'w')
	tab_string = ""
	for struct in list(database.keys()):
		struct_filename = locations['mainpath'] + locations['rpdb'] + struct + '.pdb'
		struct_file = open(struct_filename, 'r')
		text = struct_file.read().split("\n")
		struct_file.close()

		from_pdb_dict = {}
		tech_list = ['NMR', 'X-RAY', 'THEORETICAL', 'ELECTRON']
		res_ids = {}
		b_factor = {}
		b_norm = {}
		r_factor = 'NULL'
		for line in text:
			if line[0:6] == 'EXPDTA':
				for tech in tech_list:
					if tech in line:
						from_pdb_dict['TECHNIQUE'] = tech
				if 'TECHNIQUE' not in from_pdb_dict:
					from_pdb_dict['TECHNIQUE'] = 'OTHER'
			elif line[0:10] == 'REMARK   2' and 'RESOLUTION' in line:
				fields = line.split()
				for nf in range(len(fields)):
					if 'ANGSTROM' in fields[nf]:
						if fields[nf-1] != "NULL":
							from_pdb_dict['RESOLUTION'] = float(fields[nf-1])
						else:
							from_pdb_dict['RESOLUTION'] = 0
				if 'RESOLUTION' not in from_pdb_dict and from_pdb_dict['TECHNIQUE'] != 'THEORETICAL' and from_pdb_dict['TECHNIQUE'] != 'NMR':
					raise NameError("ERROR: Resolution annotation of pdb {0} is badly formatted: {1}".format(struct_filename, line))
			elif line[0:10] == 'REMARK   3' and 'FREE R VALUE' in line and 'ERROR' not in line and 'SET' not in line and (line.split()[3] == 'FREE' or line.split()[3] == 'BIN'):
				try:
					r_factor = float(line.split()[-1])
				except ValueError:
					r_factor = 'NULL'
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
				title = ''.join(line.split()[1:])

		from_pdb_dict['CHAIN'] = {}
		for chain in res_ids:
			exclude_chain = []
			print(struct, chain)
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
			if not filters['NMR'] and from_pdb_dict['TECHNIQUE'] == 'NMR':
				exclude_chain.append('NMR structure')

			# Theoretical model check
			if not filters['THM'] and from_pdb_dict['TECHNIQUE'] == 'THEORETICAL':
				exclude_chain.append('Theoretical model')

			# Holes greater than hole threshold check
			if filters['hole_thr'] > 0:
				for n in range(1,len(res_ids[chain])):
					if res_ids[chain][n][0] - res_ids[chain][n-1][0] > filters['hole_thr']:
						exclude_chain.append('Contains hole longer than {0} residues'.format(filters['hole_thr']))
						break

			# Strange residues
			for res in res_ids[chain]:
				if from3to1(res[1]) == '0':
					exclude_chain.append('One or more residues with unsupported names')
					break

			# Resolution check
			if 'resolution' in filters and from_pdb_dict['TECHNIQUE'] != 'THEORETICAL' and from_pdb_dict['TECHNIQUE'] != 'NMR':
				if from_pdb_dict['RESOLUTION'] > filters['resolution']:
					exclude_chain.append('Resolution is higher than {0}'.format(filters['resolution']))
				elif from_pdb_dict['RESOLUTION'] == 0:
					exclude_chain.append('No resolution information found'.format(filters['resolution']))


			# Check consistency with OPM
			if opm_data and struct in opm_data and chain in opm_data[struct] and chain not in database[struct][1]['CHAIN']:
				n_opm = int(opm_data[struct][chain])
				print(n_pdbtm, n_opm)
				if n_opm != n_pdbtm:
					exclude_chain.append('Chain with {0} TM domains has different number of TM domains in OPM: {1}'.format(n_pdbtm, n_opm))

			# Was there something wrong?
			if not exclude_chain:
				s_type = database[struct][1]['CHAIN'][chain][0]['TYPE']
				if struct not in instructions:
					instructions[struct] = {}
				instructions[struct][chain] = (s_type, n_pdbtm)
				chain_filename = locations['mainpath'] + locations['cpdb'] + struct + '_' + chain + '.pdb'
				chain_file = open(chain_filename, 'w')
				for line in text:
					if line[0:4] == 'ATOM':
						chain_file.write(line + '\n')
				chain_file.close()
			else:
#				print(exclude_chain)
				exclusions_file.write(struct + '_' + chain + '\t\t' + exclude_chain[0] + '\n')
				for nl in range(1, len(exclude_chain)):
					exclusions_file.write(' '*len(struct) + ' ' + ' '*len(chain) + '\t\t' + exclude_chain[nl] + '\n')

			if b_factor[chain] == 0:
				b_factor[chain] = 'NULL'
			else:
				b_factor[chain] = b_factor[chain] / b_norm[chain]
			from_pdb_dict['TITLE'] = title
			from_pdb_dict['RFACTOR'] = r_factor
			from_pdb_dict['CHAIN'][chain] = [{'CHAINID' : chain,
                                                          'AVG_BFACTOR' : b_factor[chain],
			                                  'NRES' : len(res_ids[chain])},
			                                 {'RESIDS' : [x[0] for x in res_ids[chain]],
			                                  'RESNAMES' : [from3to1(x[1]) for x in res_ids[chain]]}]

		# Introduce from_pdb_dict as a key of the database
		database[struct][1]['FROM_PDB'] = from_pdb_dict

		if not tab_string:
			for key in [x for x in sorted(list(database[struct][1]['FROM_PDB'].keys())) if x != 'CHAIN']:
				tab_string += "{0:16} ".format(key)
			tab_string += "{0:5} ".format('CHAIN')
			for key in [x for x in sorted(list(database[struct][1]['FROM_PDB']['CHAIN'][chain][0].keys())) if x != 'CHAINID']:
				tab_string += "{0:16} ".format(key)
			tab_file.write(tab_string[:-1] + "\n")

		tab_string = ""
		for key in [x for  x in sorted(list(database[struct][1]['FROM_PDB'].keys())) if x != 'CHAIN']:
			tab_string += "{0:16} ".format(database[struct][1]['FROM_PDB'][key])
		counter = 0
		strlen = len(tab_string)
		for chain in sorted(list(database[struct][1]['FROM_PDB']['CHAIN'])):
			if counter > 0:
				tab_string = " "*strlen
			tab_string += "{0:5} ".format(database[struct][1]['FROM_PDB']['CHAIN'][chain][0]['CHAINID'])
			for key in [x for x in sorted(list(database[struct][1]['FROM_PDB']['CHAIN'][chain][0].keys())) if x != 'CHAINID']:
				tab_string += "{0:16} "
			tab_file.write(tab_string[:-1] + "\n")
			counter += 1
		tab_string = "XXX"

	exclusions_file.close()
	tab_file.close()

	return database, instructions


def structure_sorter(locations, instructions):
	ssd = {'alpha' : 'a', 'beta' : 'b'}
	nprogr = {}
	for struct in instructions:
		for chain in instructions[struct]:
			ss = instructions[struct][chain][0]
			ntm = instructions[struct][chain][1]
			destination_dir = locations['mainpath'] + ss + '/' + str(ntm) + '/'
			if not os.path.exists(destination_dir):
				os.mkdir(destination_dir)
				os.mkdir(destination_dir + 'structures/')
				nprogr[ntm] = 1
			else:
				nprogr[ntm] += 1
			shutil.copy(locations['mainpath'] + locations['cpdb'] + struct + '_' + chain + '.pdb', destination_dir + 'structures/')
			code = ssd[ss] + str(ntm).zfill(3) + str(nprogr[ntm]).zfill(6)
			codes_file = open(destination_dir + 'struct_codes.dat', 'a+')
			if struct not in codes_file.read().split('\n'):
				codes_file.write(code + '\t\t' + struct)
			codes_file.close()


# Library function
def generate_chain_pdb_files(locations, database, filters):
	# Hardcoded variables
	this_name = 'genclib'
	indent = " "*len(header(this_name))
	version = 3.1

	# Checks
	for path_name in [locations['mainpath'] + locations[x] for x in list(locations.keys()) if x != 'installpath' and x != 'mainpath' and x!= 'main']:
		if not os.path.exists(path_name):
			logmsg = header(this_name) + "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name)
			write_log(this_name, logmsg)	
			raise NameError(logmsg)

	# Structure checker
	database, fs_ordering = checker(locations, database, filters)

	# Filesystem sorting
	structure_sorter(locations, fs_ordering)

	return database
