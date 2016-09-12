# Name: straln.py
# Language: python3
# Libraries: 
# Description: Runs the structural alignments needed for HOMEP
# Author: Edoardo Sarti
# Date: Aug 17 2016

import os, sys, multiprocessing, subprocess, re, time

def repo_inspector(repo_filename):
	repo_file = open(repo_filename, 'r')
	text = repo_file.read().split('\n')
	repo_file.close()
	repo_info = []
	for line in text:
		fields = line.split()
		if line and fields[0] == 'BEGIN':
			record = True
			recorded_text = ""
		if not record and not line:
			continue
		if record:
			recorded_text = line + '\n'
		if fields[0] == 'CHAIN_1:':
			chain_1 = fields[1]
		if fields[0] == 'CHAIN_2:':
			chain_2 = fields[1]
		if fields[0] == 'END':
			if chain_1 not in repo_info:
				repo_info[chain_1] = {}
			if chain_2 in repo_info[chain_1]:
				print("WARNING: multiple occurrences of couple {0} {1} in file {2}".format(chain_1, chain_2, repo_filename))
			repo_info[chain_1][chain_2] = recorded_text
			record = False
	return repo_info


def FrTMjob(data):
	locations, target, exelist = data
	clade, topology, chain_1, straln_path = target
	
	topology_path = locations['FSYS']['mainpath'] + clade + '/' + topology + '/'
	sequence_path = topology_path + locations['FSYS']['TREE']['seqaln'] + 'seq_' + chain_1 + '.dat'
	structure_path = topology_path + locations['FSYS']['TREE']['straln'] + 'str_' + chain_1 + '.dat'
	seq_repo_path = locations['FSYS']['mainpath'] + locations['FSYS']['repochains'] + locations['FSYS']['TREE']['seqaln'] 
	str_repo_path = locations['FSYS']['mainpath'] + locations['FSYS']['repochains'] + locations['FSYS']['TREE']['straln']

	# Checks for needed locations:
	# Superfamily path
	if not os.path.exists(topology_path):
		raise NameError("ERROR: Superfamily {0} not found in path {1}".format(clade+' '+topology, topology_path)
	# structure/ folder
	if not os.path.exists(topology_path + locations['FSYS']['TREE']['str']):
		raise NameError("ERROR: {0} folder not found in path {1}".format(locations['FSYS']['TREE']['str'], topology_path))
	# Main pdb file
	if not os.path.exists(topology_path + locations['FSYS']['TREE']['str'] + chain_1 + '.pdb'):
		raise NameError("ERROR: File {0} not found in {1}".format(chain_1 + '.pdb', topology_path + locations['FSYS']['TREE']['str']))
	# Secondary pdb files
	for chain_2 in exelist:
		if not os.path.exists(topology_path + locations['FSYS']['TREE']['str'] + chain_2 + '.pdb'):
			raise NameError("ERROR: File {0} not found in {1}".format(chain_2 + '.pdb', topology_path + locations['FSYS']['TREE']['str']))
	
	# Creates, if needed, the alignment locations
	for n, x in enumerate(locations['FSYS']['TREE'].items()):
		if n > 0 and not os.path.exists(topology_path + x[1]):
			os.mkdir(topology_path + x[1])

	# Creates, if needed, the repository locations
	if not os.path.exists(seq_repo_path):
		os.mkdir(seq_repo_path)
	if not os.path.exists(str_repo_path):
		os.mkdir(str_repo_path)

	# Checks if the main sequence file already exists. If it does, checks again if none of the alignments in exelist
	# is contained in the file. Cancels from exelist any found alignment.
	if os.path.exists(sequence_path):
		sequence_file = open(sequence_path, 'r')
		text = sequence_file.read().split('\n')
		for line in text:
			if 'CHAIN_2:' in line and field[1] in exelist:
				logmsg = "Update repository with {0}_{1} sequence alignment".format(exelist[2], field[1])
				update_repository.append(field[1])
				exelist.remove[field[1]]
		sequence_file.close()

	# Checks if the main structure file already exists. If it does, checks again if none of the alignments in exelist
	# is contained in the file. Cancels from exelist any found alignment.
	if os.path.exists(structure_path):
		structure_file = open(structure_path, 'r')
		text = structure_file.read().split('\n')
		for line in text:
			if 'CHAIN_2:' in line and field[1] in exelist:
				logmsg = "Update repository with {0}_{1} alignment".format(exelist[2], field[1])
				update_repository.append(field[1])
				exelist.remove[field[1]]
		structure_file.close()
	
	# Creates the temporary folder for sequence alignments
	aln_tmpfolder_path = topology_path + 'alignments/fasta/tmp_' + chain_1 + '/'
	if not os.path.exists(aln_tmpfolder_path):
		os.mkdir(aln_tmpfolder_path)

	# Creates the temporary folder for structure alignments
	straln_tmpfolder_path = topology_path + 'alignments/str_alns/tmp_' + chain_1 + '/'
	if not os.path.exists(straln_tmpfolder_path):
		os.mkdir(straln_tmpfolder_path)

	pdb1_filename = topology_path + 'structures/' + chain_1 + '.pdb'
	for chain_2 in exelist:
		# Defines filenames
		pdb2_filename = topology_path + 'structures/' + chain_2 + '.pdb'
		seq_output_filename = aln_tmpfolder_path + 'aln_' + chain_1 + '_' + chain_2 + '.tmp'
		str_output_filename = 'straln_' + chain_1 + '_' + chain_2 + '.tmp' # This filename is without path for a constraint on flags length of frtmalign (see below)
		stdout_filename = aln_tmpfolder_path + 'output_' + chain_1 + '_' + chain_2 + '.tmp'

		# If sequence and structure temporary files are present, jump
		if os.path.exists(seq_output_filename) and os.path.exists(str_output_filename):
			continue

		# Fr-TM-align call
		# The Fr-TM-align command cannot exceed a certain length. Thus, we are compelled to run it locally into structures/
		# The stdout file from Fr-TM-align can go directly to the right temporary folder (but needs to be cleaned)
		stdout_file = open(stdout_filename, 'w')
		fnull = open(os.devnull, 'w')
		print(straln_path, pdb1_filename[-10:], pdb2_filename[-10:], '-o', str_output_filename)
		p = subprocess.Popen([straln_path, pdb1_filename[-10:], pdb2_filename[-10:], '-o', str_output_filename], stdout=stdout_file, stderr=fnull, cwd=topology_path+'structures/')
		p.wait()
		fnull.close()
		stdout_file.close()

		# Moves the Fr-TM-align output file into the structure temporary folder
		os.rename(topology_path + 'structures/' + str_output_filename, straln_tmpfolder_path + str_output_filename)

		# Reads and then removes the stdout file from Fr-TM-align
		stdout_file = open(stdout_filename, 'r')
		text = stdout_file.read().split('\n')
		stdout_file.close()
		os.remove(stdout_filename)
		
		# From the stdout file from Fr-TM-align, it takes the RMSD, TM-score and the two aligned sequences
		chkaln = -1000
		for nl in range(len(text)):
			if "Aligned length" in text[nl]:
				fields = re.split('=|,|\s',text[nl])
				fields = list(filter(None, fields))
				RMSD = float(fields[4])
				TMscore = float(fields[6])
			elif chkaln+1 == nl:
				seq_1 = text[nl]
			elif chkaln+3 == nl:
				seq_2 = text[nl]
			elif "denotes the residue pairs of distance" in text[nl]:
				chkaln = nl
		
		# Creates a sequence temporary file already correctly formatted
		seq_output_file = open(seq_output_filename, 'w')
		seq_output_file.write(">" + chain_1 + "\n" + seq_1.replace('\x00', '') + "\n>" + chain_2 + "\n" + seq_2.replace('\x00', '') + "\n\nRMSD\t{0:.2f}\nTM-score\t{1:.5f}\n\n".format(RMSD, TMscore))
		seq_output_file.close()

	# Writes on the main sequence file. Each alignment begins with a "BEGIN" line, followed by two lines reporting the two chain names
	# (format: "CHAIN_X: chain_name", where X={1,2}), and ends with an "END" line.
	sequence_file = open(sequence_path, 'a')
	for tmp_filename in sorted(os.listdir(aln_tmpfolder_path)):
		sequence_file.write("BEGIN \nCHAIN_1: " + chain_1  + "\nCHAIN_2: " + tmp_filename[-10:-4] + "\n")
		tmp_file = open(aln_tmpfolder_path + tmp_filename)
		text = tmp_file.read().split('\n')
		tmp_file.close()
		for line in text:
			sequence_file.write(line+'\n')
		sequence_file.write("END\n\n\n")
		os.remove(aln_tmpfolder_path + tmp_filename)
	time.sleep(1)
	os.rmdir(aln_tmpfolder_path)
	sequence_file.close()

#	repo_inspector()

	# Writes on the main structure file. Each alignment begins with a "BEGIN" line, followed by two lines reporting the two chain names
	# (format: "CHAIN_X: chain_name", where X={1,2}), and ends with an "END" line.
	structure_file = open(structure_path, 'a')
	for tmp_filename in sorted(os.listdir(straln_tmpfolder_path)):
		structure_file.write("BEGIN \nCHAIN_1: " + chain_1  + "\nCHAIN_2: " + tmp_filename[-10:-4] + "\n")
		tmp_file = open(straln_tmpfolder_path + tmp_filename)
		text = tmp_file.read().split('\n')
		tmp_file.close()
		for line in text:
			structure_file.write(line+'\n')
		structure_file.write("END\n\n\n")
		os.remove(straln_tmpfolder_path + tmp_filename)
	time.sleep(1)
	os.rmdir(straln_tmpfolder_path)
	structure_file.close()
	

def calculate_seqid(alignment):
	ntot = 0
	naln = 0
	for na in range(len(alignment[0])):
		if alignment[0][na] != '-' and alignment[1][na] != '-':
			ntot += 1
		if alignment[0][na] == alignment[1][na]:
			naln += 1
	if ntot > 0:
		return naln/ntot
	else:
		return 0


def make_new_table(locations, external_filename):
	names = {'alpha' : 'a', 'beta' : 'b'}

	topologies = []
	for ss in 'alpha', 'beta':
		for i in os.listdir(locations['FSYS']['mainpath'] + ss + '/'):
			if re.match('^\d*$', str(i)):
				topologies.append((ss, i))

	table_filename = locations['FSYS']['mainpath'] + external_filename
	table_file = open(table_filename, 'w')
	table = {}
	for sf in sorted(topologies, key = lambda x: (x[0], int(x[1]))):
		if not sf[0] in table:
			table[sf[0]] = {}
		table[sf[0]][sf[1]] = {}
		topology_seqaln_path = locations['FSYS']['mainpath'] + sf[0] + '/' + sf[1] + '/alignments/fasta/'
		if os.path.exists(topology_seqaln_path): 
			files_in_seqaln_path = os.listdir(topology_seqaln_path)
		else:
			files_in_seqaln_path = []
		for seqaln_filename in files_in_seqaln_path:
			if seqaln_filename[0:4] == 'seq_' and seqaln_filename[-4:] == '.dat':
				seqaln_file = open(topology_seqaln_path + seqaln_filename, 'r')
				text = seqaln_file.read().split('\n')
				seqaln_file.close()
				for nline in range(len(text)):
					if not text[nline]:
						continue
					fields = text[nline].split()
#					print(seqaln_filename, fields)
					if fields[0] == 'BEGIN':
						continue
					elif fields[0] == 'CHAIN_1:':
						chain_1 = fields[1]
					elif fields[0] == 'CHAIN_2:':
						chain_2 = fields[1]
					elif fields[0] == '>' + chain_1:
						seq_1 = text[nline+1]
					elif fields[0] == '>' + chain_2:
						seq_2 = text[nline+1]
					elif fields[0] == 'RMSD':
						RMSD = float(fields[1])
					elif fields[0] == 'TM-score':
						tmscore = float(fields[1])
					elif fields[0] == 'END':
						if not (chain_1 and chain_2 and seq_1 and seq_2):
							raise NameError("ERROR: file is corrupted")
						seqid = calculate_seqid((seq_1, seq_2))
						table_file.write("{0}\t{1}\t{2}\t{3}\t{4:10.8f}\t{5:10.8f}\t{6:10.6f}\t\t{7}\n".format(names[sf[0]], str(int(sf[1])).zfill(3), chain_1, chain_2, seqid, tmscore, RMSD, topology_seqaln_path+seqaln_filename))
						if chain_1 not in table[sf[0]][sf[1]]:
							table[sf[0]][sf[1]][chain_1] = {}
						table[sf[0]][sf[1]][chain_1][chain_2] = (seqid, tmscore, RMSD, topology_seqaln_path+seqaln_filename)
	table_file.close()
	return table
							

def structure_alignment(options, locations):
	aligner_path = options['straln_path']
	np = int(options['number_of_procs'])
	external_filename = options['output_tab']

	already_processed = []
	ex_list = {}
	hidden_repository_filename = locations['FSYS']['mainpath'] + '.structure_alignments.dat'
	if os.path.exists(hidden_repository_filename):
		hidden_repository_file = open(hidden_repository_filename, 'r')
		text = hidden_repository_file.read().split("\n")
		for line in text:
			if not line:
				continue
			fields = line.split()
			already_processed.append((fields[2], fields[3]))
		hidden_repository_file.close()
	repository_filename = locations['FSYS']['mainpath'] + external_filename
	if os.path.exists(repository_filename):
		repository_file = open(repository_filename, 'r')
		text = repository_file.read().split("\n")
		for line in text:
			if not line:
				continue
			fields = line.split()
			if not (fields[2], fields[3]) in already_processed:
				already_processed.append((fields[2], fields[3]))
		external_filename = 'new_' + external_filename
		repository_file.close()
	
	topologies = []
	for ss in 'alpha', 'beta':
		for i in os.listdir(locations['FSYS']['mainpath'] + ss + '/'):
			if re.match('^\d*$', str(i)):
				topologies.append((ss, i))

	ex_check = {}
	exelist_filename = locations['FSYS']['mainpath'] + '.scheduled_alignments.dat'
	if os.path.exists(exelist_filename):
		exelist_file = open(exelist_filename, 'r')
		text = exelist_file.read().split('\n')
		for line in text:
			if not line:
				continue
			fields = line.split()
			if fields[2] not in ex_list:
				ex_list[fields[2]] = []
			ex_list[fields[2]].append(fields[3])
		exelist_file.close()

	exelist = {}
	exelist_file = open(exelist_filename, 'a')
	for sf in topologies:
		structs = [x[-10:-4] for x in os.listdir(locations['FSYS']['mainpath'] + sf[0] + '/' + sf[1] + '/structures/') if x[-4:] == '.pdb']
#		print(structs)
		if len(structs) > 1:
			for s1 in structs:
				exelist[s1] = {}
				for s2 in structs:
					if s1 == s2 or (s1, s2) in already_processed:
						continue
#					print("sf0", sf[0], "sf1", sf[1], "s1", s1, "s2", s2)
					exelist[s1][s2] = (sf[0], sf[1], s1, s2)
					if s1 not in ex_list or s2 not in ex_list[s1]:
						exelist_file.write("{0}\t{1}\t{2}\t{3}\n".format(sf[0], sf[1], s1, s2))
	exelist_file.close()

#	print("exelist", sorted(list(exelist.keys())))	

	data = []
	for i in sorted(list(exelist.keys())):
		exesublist = []
		if not list(exelist[i].keys()):
			continue
		for j in sorted(list(exelist[i].keys())):
			exesublist.append(exelist[i][j][3])
		jtmp = sorted(list(exelist[i].keys()))[0]
		data.append((locations, (exelist[i][jtmp][0], exelist[i][jtmp][1], exelist[i][jtmp][2], aligner_path), exesublist))

	pool = multiprocessing.Pool(processes=np)
	pool_outputs = pool.map(FrTMjob, data)
	pool.close()
	pool.join()

	table = make_new_table(locations, external_filename)

	return table


def check_runs():
	instructions_filename = locations['FSYS']['mainpath'] + '.topology_classification.txt'
	instructions_file = open(instructions_filename, 'r')
	text = instructions_file.read().split('\n')
	instructions_file.close()


"""
import genfsys, genrlib, genclib, clusterize
main_dir = '/u/esarti/LoBoS_Workspace_ES/projects/HOMEP/HOMEP3.1/scripts/generate_library_PDBTM/test/HOMEP_3.1_2016_09_01/'
straln_path = '/v/apps/csb/frtmalign/frtmalign.exe'
pdbtm_file_path = 'pdbtmall_reduced'
output_tab = 'structure_alignments.dat' 
HOMEP_filename = 'HOMEP_prova.dat'
seqid_thr = 0.85
tmscore_thr = 0.6
np = 8
filters = {'resolution' : 3.5,
           'NMR' : False,
           'THM' : False,
           'hole_thr' : 100}


locations = genfsys.filesystem_info(str(main_dir))
pdbtm_data = genrlib.generate_raw_pdb_library(locations, pdbtm_file_path)
pdbtm_data = genclib.generate_chain_pdb_files(locations, pdbtm_data, filters)
print("init straln")
table = structure_alignment(locations, straln_path, np, str(output_tab))
o = clusterize.clusterize(locations, pdbtm_data, table, '', HOMEP_filename, seqid_thr, tmscore_thr)
#print(table)
"""
