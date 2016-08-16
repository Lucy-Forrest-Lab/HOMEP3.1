# Name: genclib.py
# Language: python3
# Libraries:
# Description: Generates HOMEP chain library
# Author: Edoardo Sarti
# Date: Aug 15 2016

import ppm_segments 

def from3to1(resname):
	3to1 = {'ALA' : 'A',
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

	if resname in list(3to1.keys()):
		return 3to1[resname]
	else:
		return '0'


# Structure checker
def checker(locations, database, filters):
	opm_data = {}
	if 'OPM_TMdoms' in filters:
		opm_data = ppm_segments.OPM_TMdoms()

	instructions = {}
	for struct in list(database.keys()):
		struct_filename = locations['raw_pdbs'] + struct + '.pdb'
		struct_file = open(struct_filename, 'r')
		text = struct_file.read().split("\n")
		struct_file.closed()

		from_pdb_dict = {}
		tech_list = ['NMR', 'X-RAY', 'THEORETICAL', 'ELECTRON']
		res_ids = []
		b_factor = {}
		b_norm = {}
		for line in text:
			if line[0:6] == 'EXPDTA':
				for tech in tech_list:
					if tech in line:
						from_pdb_dict['TECHNIQUE'] = tech
				if 'TECHNIQUE' not in from_pdb_dict:
					from_pdb_dict['TECHNIQUE'] = 'OTHER'
			if line[0:10] == 'REMARK   2' and 'RESOLUTION' in line:
				fields = line.split()
				for nf in range(len(fields)):
					if 'ANGSTROM' in fields[nf]:
						from_pdb_dict['RESOLUTION'] = float(fields[nf-1])
				if 'RESOLUTION' not in from_pdb_dict:
					raise NameError("ERROR: Resolution annotation of pdb {0} is badly formatted: {1}".format(struct_filename, line))
			if line[0:4] == 'ATOM':
				if not line[21]:
					raise NameError("ERROR: There is an ATOM without chain name: {0}".format(line))
				if line[21] not in res_ids:
					res_ids[line[21]] = []
					b_factor[line[21]] = 0
					b_norm[line[21]] = 0
				if res_ids[line[21]] and res_ids[line[21]][0] != int(line[22:26]):
#					if int(line[22:26]) - res_ids[line[21]][0] > filters['max_hole_length']:
#						exclude = True
					res_ids[line[21]].append((int(line[22:26]), line[17:20]))
				if line[60:66]:
					b_factor[line[21]] += float(line[60:66])
					b_norm[line[21]] += 1

		from_pdb_dict['CHAIN'] = {}
		exclusions_filename = locations['pdbs'] + 'exclusions.txt'
		exclusions_file(exclusions_filename, 'w')
		for chain in res_ids:
			exclude_chain = []
			n_pdbtm = database[struct][1]['CHAIN'][chain][0]['NUM_TM']
			
			# NMR check
			if not 'NMR' in filters and from_pdb_dict['TECHNIQUE'] == 'NMR':
				exclude_chain.append('NMR structure')

			# Holes greater than hole threshold check
			if filters['max_hole_length'] > 0:
				for n in range(1,len(res_ids[chain])):
					if res_ids[chain][n] - res_ids[chain][n-1] > filters['max_hole_length']:
						exclude_chain.append('Contains hole longer than {0} residues'.format(filters['max_hole_length']))
						break

			# Strange residues
			for res in res_ids[chain]:
				if from3to1(res[1]) == '0':
					exclude_chain.append('One or more residues with unsupported names')
					break

			# Resolution check
			if 'resolution' in filters:
				if from_pdb_dict['RESOLUTION'] > filters['resolution']:
					exclude_chain.append('Resolution is higher than {0}'.format(filters['resolution']))


			# Check consistency with OPM
			if opm_data and opm_data[struct] and opm_data[struct][chain]:
				n_opm = opm_data[struct][chain]
				if n_opm != n_pdbtm:
					exclude_chain.append('Chain with {0} TM domains has different number of TM domains in OPM: {1}'.format(n_pdbtm, n_opm)

			# Was there something wrong?
			if not exclude_chain:
				s_type = database[struct][0]['TYPE']
				if struct not in instructions:
					instructions[struct] = {}
				instructions[struct][chain] = (s_type, n_pdbtm)
				chain_filename = locations['pdbs'] + struct + '_' + chain + '.pdb'
				chain_file = open(chain_filename, 'w')
				for line in text:
					if line[0:4] == 'ATOM':
						chain_file.write(line + '\n')
				chain_file.close()
			else:
				exclusions_file.write(struct + '_' + chain + '\t\t' + exclude_chain + '\n')
				for nl in range(1, len(exclude_chain)):
					exclusions_file.write(' '*len(struct) + ' ' + ' '*len(chain) + '\t\t' + exclude_chain[nl] + '\n')

			from_pdb_dict['CHAIN'][chain] = [{'CHAINID' : chain,
                                                        'AVG_BFACTOR' : b_factor[chain] / b_norm[chain],
			                                'NRES' : len(res_ids[chain])},
			                               {'RESIDS' : [x[0] for x in res_ids[chain]],
			                                'RESNAMES' : [from3to1(x[1]) for x in res_ids[chain]]}]

		exclusions_file.close()

		# Introduce from_pdb_dict as a key of the database
		database[struct][1]['FROM_PDB'] = from_pdb_dict

	return database, instructions


def structure_sorter(locations, fs_ordering):
	ssd = {'alpha' : 'a', 'beta' : 'b'}
	nprogr = {}
	for struct in instructions:
		for chain in instructions[struct]:
			ss = fs_ordering[struct][chain][0]
			ntm = fs_ordering[struct][chain][1]
			destination_dir = locations['main'] + ss + '/' + ntm
			if not os.path.esxists(destination_dir):
				os.makedir(destination_dir)
				os.makedir(destination_dir + '/structures')
				nprogr[ntm] = 1
			else:
				nprogr[ntm] += 1
			shutil.copy(locations['pdbs'] + struct + '_' + chain + '.pdb', destination_dir + '/structures')
			code = ssd[ss] + str(ntm).zfill(3) + str(nprogr[ntm]).zfill(6)
			codes_file = open(destination_dir + 'struct_codes.dat', 'a')
			codes_file.write(code + '\t\t' + struct)
			codes_file.close()


# Library function
def generate_chain_pdb_files(locations, database, filters):
	# Hardcoded variables
	this_name = 'genclib'
	indent = " "*len(header(this_name))
	version = 3.1

	# Checks
	for path_name in [locations['installpath']+x for x in locations if x != 'installpath']:
		if not os.path.exists(path_name):
			logmsg = header(this_name) + "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name)
			write_log(this_name, logmsg)	
			raise NameError(logmsg)

	# Structure checker
	database, fs_ordering = checker(locations, database, filters)

	# Filesystem sorting
	structure_sorter(locations, fs_ordering)

	return database
