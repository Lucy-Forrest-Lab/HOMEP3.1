# Name: straln.py
# Language: python3
# Libraries: 
# Description: Runs the structural alignments needed for HOMEP
# Author: Edoardo Sarti
# Date: Aug 17 2016

import os
import sys
import multiprocessing
import subprocess
import re
import time
import Bio

def repo_inspector(repo_filename):
	repo_file = open(repo_filename, 'r')
	text = repo_file.read().split('\n')
	repo_file.close()
	repo_info = {}
	for line in text:
		fields = line.split()
		if line and fields[0] == 'BEGIN':
			chain_1 = fields[2]
			chain_2 = fields[4]
			record = True
		if not record and not line:
			continue
		if record:
			recorded_text = line + '\n'
		if fields[0] == 'END':
			if chain_1 not in repo_info:
				repo_info[chain_1] = {}
			if chain_2 in repo_info[chain_1]:
				print("WARNING: multiple occurrences of couple {0} {1} in file {2}".format(chain_1, chain_2, repo_filename))
			repo_info[chain_1][chain_2] = recorded_text
			record = False
	return repo_info


def write_on_repo(repo_filename, textdict, append=False):
	if append:
		repo_file = open(repo_filename, 'a')
	else:
		repo_file = open(repo_filename, 'w')
	for val1 in textdict:
		for val2 in textdict[val1]:
			repo_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n" + textdict + "\nEND\n\n\n")
	repo_file.close()


def FrTMjob(data):
	locations, target, exelist = data
	topologytype, topology, chain_1, straln_path = target
	
	topology_path = locations['FSYSPATH'][topologytype] + topology + '/'
	strfasta_filename = topology_path + locations['TREE']['straln'] + 'str_' + chain_1 + '_fasta.dat'
	strpdb_filename = topology_path + locations['TREE']['straln'] + 'str_' + chain_1 + '_pdb.dat'
	str_repo_path = locations['FSYSPATH']['repocstraln']
	seq_repo_path = locations['FSYSPATH']['repocseqaln']

	# Checks for needed locations:
	# Superfamily path
	if not os.path.exists(topology_path):
		raise NameError("ERROR: Superfamily {0} not found in path {1}".format(topologytype+' '+topology, topology_path))
	# structure/ folder
	if not os.path.exists(topology_path + locations['TREE']['str']):
		raise NameError("ERROR: {0} folder not found in path {1}".format(locations['TREE']['str'], topology_path))
	# Main pdb file
	if not os.path.exists(topology_path + locations['TREE']['str'] + chain_1 + '.pdb'):
		raise NameError("ERROR: File {0} not found in {1}".format(chain_1 + '.pdb', topology_path + locations['TREE']['str']))
	# Secondary pdb files
	for chain_2 in exelist:
		if not os.path.exists(topology_path + locations['TREE']['str'] + chain_2 + '.pdb'):
			raise NameError("ERROR: File {0} not found in {1}".format(chain_2 + '.pdb', topology_path + locations['TREE']['str']))
	
	# Creates, if needed, the alignment locations
	for name, val in locations['TREE'].items():
		if 'aln' in name and not os.path.exists(topology_path + val):
			os.mkdir(topology_path + val)

	# Creates, if needed, the repository locations
	if not os.path.exists(seq_repo_path):
		os.mkdir(seq_repo_path)
	if not os.path.exists(str_repo_path):
		os.mkdir(str_repo_path)

	repo_info_str_pdb = repo_inspector(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_pdb.dat')
	repo_info_str_fasta = repo_inspector(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_fasta.dat')
	repo_info_seq_fasta = repo_inspector(locations['FSYSPATH']['repocseqaln'] + 'seq_' + chain_1 + '_fasta.dat')

	"""
	# Checks if the main sequence file already exists. If it does, checks again if none of the alignments in exelist
	# is contained in the file. Cancels from exelist any found alignment.
	if os.path.exists(fasta_path):
		fasta_file = open(fasta_path, 'r')
		text = fasta_file.read().split('\n')
		for line in text:
			if 'CHAIN_2:' in line and field[1] in exelist:
				logmsg = "Update repository with {0}_{1} sequence alignment".format(exelist[2], field[1])
				update_repository.append(field[1])
				exelist.remove[field[1]]
		fasta_file.close()

	# Checks if the main structure file already exists. If it does, checks again if none of the alignments in exelist
	# is contained in the file. Cancels from exelist any found alignment.
	if os.path.exists(pdb_path):
		pdb_file = open(pdb_path, 'r')
		text = pdb_file.read().split('\n')
		for line in text:
			if 'CHAIN_2:' in line and field[1] in exelist:
				logmsg = "Update repository with {0}_{1} alignment".format(exelist[2], field[1])
				update_repository.append(field[1])
				exelist.remove[field[1]]
		pdb_file.close()
	"""

	# Creates the temporary folder for sequence alignments
	fasta_tmpfolder_path = topology_path + locations['TREE']['stralns'] + 'tmp_' + chain_1 + '_fasta/'
	if not os.path.exists(fasta_tmpfolder_path):
		os.mkdir(fasta_tmpfolder_path)

	# Creates the temporary folder for structure alignments
	pdb_tmpfolder_path = topology_path + locations['TREE']['stralns'] + 'tmp_' + chain_1 + '_pdb/'
	if not os.path.exists(pdb_tmpfolder_path):
		os.mkdir(pdb_tmpfolder_path)

	pdb1_filename = topology_path + locations['TREE']['str'] + chain_1 + '.pdb'
	for chain_2 in exelist:
		if repo_info_str_fasta[chain_1] and repo_info_str_fasta[chain_1][chain_2] and repo_info_str_pdb[chain_1] and repo_info_str_pdb[chain_1][chain_2]:
			continue

		# Defines filenames
		pdb2_filename = topology_path + locations['TREE']['str'] + chain_2 + '.pdb'
		fasta_output_filename = fasta_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_fasta.tmp'
		pdb_output_filename = 'straln_' + chain_1 + '_' + chain_2 + '_pdb.tmp' # This filename is without path for a constraint on flags length of frtmalign (see below)
		stdout_filename = fasta_tmpfolder_path + 'output_' + chain_1 + '_' + chain_2 + '.tmp'

		# If sequence and structure temporary files are present, jump
		if os.path.exists(fasta_output_filename) and os.path.exists(pdb_output_filename):
			continue

		# Fr-TM-align call
		# The Fr-TM-align command cannot exceed a certain length. Thus, we are compelled to run it locally into structures/
		# The stdout file from Fr-TM-align can go directly to the right temporary folder (but needs to be cleaned)
		stdout_file = open(stdout_filename, 'w')
		fnull = open(os.devnull, 'w')
		print(straln_path, pdb1_filename[-10:], pdb2_filename[-10:], '-o', pdb_output_filename)
		p = subprocess.Popen([straln_path, pdb1_filename[-10:], pdb2_filename[-10:], '-o', pdb_output_filename], stdout=stdout_file, stderr=fnull, cwd=topology_path+locations['TREE']['str'])
		p.wait()
		fnull.close()
		stdout_file.close()

		# Moves the Fr-TM-align output file into the structure temporary folder
		os.rename(topology_path + locations['TREE']['str'] + pdb_output_filename, pdb_tmpfolder_path + pdb_output_filename)

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
		fasta_output_file = open(fasta_output_filename, 'w')
		fasta_output_file.write(">" + chain_1 + "\n" + seq_1.replace('\x00', '') + "\n>" + chain_2 + "\n" + seq_2.replace('\x00', '') + "\n\nRMSD\t{0:.2f}\nTM-score\t{1:.5f}\nstr_SEQID\t{2:.5f}\n\n".format(RMSD, TMscore, calculate_seqid((seq_1, seq_2))))
		fasta_output_file.close()

	fasta1_filename = topology_path + locations['TREE']['seq'] + chain_1 + '.fa'
	fasta1_file = open(fasta1_filename, 'r')
	text = fasta1_file.read().split('\n')
	for line in text:
		if line and line[0] != '>':
			sequence_1 = line.strip()
			break
	fasta1_file.close()
	seqfasta_file = open(seqfasta_filename, 'w')
	for chain_2 in exelist:
		seqfasta_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n")
		if not (repo_info_seq_fasta[chain_1] and repo_info_seq_fasta[chain_1][chain_2]):
			fasta2_filename = topology_path + locations['TREE']['seq'] + chain_2 + '.fa'
			fasta2_file = open(fasta2_filename, 'r')
			text = fasta2_file.read().split('\n')
			fasta2_file.close()
			for line in text:
				if line and line[0] != '>':
					sequence_2 = line.strip()
					break
			result = Bio.pairwise2.align.globalds(sequence_1, sequence_2, Bio.SubsMat.MatrixInfo.blosum62, -10.0, -0.5)
			seqaln = [result[0][0][0], result[0][1][0]]
			text = "{0}\n{1}\n{2}\n{3}\n\nseq_SEQID: {4}\n".format(chain_1, seqaln[0], chain_2, seqaln[1], calculate_seqid(seqaln))
			seqfasta_file.write(text)
			repo_info_seq_fasta[chain_1][chain_2] = text
		else:
			for line in repo_info_seq_fasta[chain_1][chain_2]:
				seqfasta_file.write(line+"\n")
		seqfasta_file.write("END\n\n\n")
	write_on_repo(locations['FSYSPATH']['repocseqaln'] + 'seq_' + chain_1 + '_fasta.dat', repo_info_seq_fasta)

	# Writes on the main sequence file. Each alignment begins with a "BEGIN" line, followed by two lines reporting the two chain names
	# (format: "CHAIN_X: chain_name", where X={1,2}), and ends with an "END" line.
	strfasta_file = open(strfasta_filename, 'w')
	for chain_2 in exelist:
		strfasta_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n")
		if not (repo_info_str_fasta[chain_1] and repo_info_str_fasta[chain_1][chain_2]):
			tmp_filename = fasta_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_fasta.tmp'
			if not os.path.exists(tmp_filename):
				continue
			tmp_file = open(tmp_filename)
			text = tmp_file.read().split('\n')
			tmp_file.close()
			for line in text:
				strfasta_file.write(line+'\n')
			repo_info_str_fasta[chain_1][chain_2] = text
			os.remove(tmp_filename)
		else:
			for line in repo_info_str_fasta[chain_1][chain_2]:
				strfasta_file.write(line + '\n')
		strfasta_file.write("END\n\n\n")
	write_on_repo(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_fasta.dat', repo_info_str_fasta)
	time.sleep(1)
	os.rmdir(strfasta_tmpfolder_path)
	strfasta_file.close()

	# Writes on the main structure file. Each alignment begins with a "BEGIN" line, followed by two lines reporting the two chain names
	# (format: "CHAIN_X: chain_name", where X={1,2}), and ends with an "END" line.
	strpdb_file = open(strpdb_filename, 'w')
	for chain_2 in exelist:
		strpdb_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n")
		if not (repo_info_str_pdb[chain_1] and repo_info_str_pdb[chain_1][chain_2]):
			tmp_filename = pdb_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_pdb.tmp'
			if not os.path.exists(tmp_filename):
				continue
			tmp_file = open(tmp_filename)
			text = tmp_file.read().split('\n')
			tmp_file.close()
			for line in text:
				strpdb_file.write(line+'\n')
			repo_info_str_pdb[chain_1][chain_2] = text
			os.remove(tmp_filename)
		else:
			for line in repo_info_str_pdb[chain_1][chain_2]:
				strpdb_file.write(line + '\n')
		strpdb_file.write("END\n\n\n")
	write_on_repo(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_pdb.dat', repo_info_str_pdb[chain_1][chain_2])
	time.sleep(1)
	os.rmdir(pdb_tmpfolder_path)
	strpdb_file.close()

	return (repo_info_seq_fasta, repo_info_str_fasta, repo_info_str_pdb)
	

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


def make_new_table(locations, fasta_repos, external_filename):
	seq_fasta_repo, str_fasta_repo = fasta_repos
	names = {'alpha' : 'a', 'beta' : 'b'}

	instructions_filename = locations['SYSFILES']['H_topologytype']
	instructions_file = open(instructions_filename, 'r')
	text = instructions_file.read().split('\n')
	instructions_file.close()

	instructions = {}
	for line in text:
		if not line:
			continue
		fields = line.split()
		if fields[0] not in instructions:
			instructions[fields[0]] = {}
		instructions[fields[0]][fields[1]] = fields[2]

	table_filename = locations['FSYSPATH']['main'] + external_filename
	table_file = open(table_filename, 'w')
	table = {}
	for toptype in list(instructions.keys()):
		for top in instructions[toptype]:
			# REWRITE IT HERE WITH THE TWO REPO DICTIONARIES


	topologies = []
	for ss in 'alpha', 'beta':
		for i in os.listdir(locations['FSYSPATH'][ss]):
			if re.match('^\d*$', str(i)):
				topologies.append((ss, i))

	table_filename = locations['FSYSPATH']['main'] + external_filename
	table_file = open(table_filename, 'w')
	table = {}
	for sf in sorted(topologies, key = lambda x: (x[0], int(x[1]))):
		if not sf[0] in table:
			table[sf[0]] = {}
		table[sf[0]][sf[1]] = {}
		topology_strfasta_path = locations['FSYSPATH'][sf[0]] + sf[1] + '/' +  locations['TREE']['straln']
		if os.path.exists(topology_strfasta_path): 
			files_in_strfasta_path = os.listdir(topology_strfasta_path)
		else:
			files_in_strfasta_path = []
		for strfasta_filename in files_in_strfasta_path:
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
						chain_1 = fields[2]
						chain_2 = fields[4]
						continue
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
	hidden_repository_filename = locations['SYSFILES']['repocstraln']
	if os.path.exists(hidden_repository_filename):
		hidden_repository_file = open(hidden_repository_filename, 'r')
		text = hidden_repository_file.read().split("\n")
		for line in text:
			if not line:
				continue
			fields = line.split()
			already_processed.append((fields[2], fields[3]))
		hidden_repository_file.close()
	repository_filename = locations['FSYSPATH']['main'] + external_filename
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
		for i in os.listdir(locations['FSYSPATH'][ss]):
			if re.match('^\d*$', str(i)):
				topologies.append((ss, i))

	ex_check = {}
	exelist_filename = locations['SYSFILES']['H_scheduledalns']
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
		structs = [x[-10:-4] for x in os.listdir(locations['FSYSPATH'][sf[0]] + sf[1] + '/' + locations['TREE']['str']) if x[-4:] == '.pdb']
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

	repo_info_seq_fasta, repo_info_str_fasta, repo_info_str_pdb = pool_outputs[0]
	tot_repo_info_seq_fasta = repo_info_seq_fasta.copy()
	tot_repo_info_str_fasta = repo_info_str_fasta.copy()
	tot_repo_info_str_pdb = repo_info_str_pdb.copy()
	for triplet in pool_outputs[1:]:
		repo_info_seq_fasta, repo_info_str_fasta, repo_info_str_pdb = triplet
		tot_repo_info_seq_fasta.update(repo_info_seq_fasta)
		tot_repo_info_str_fasta.update(repo_info_str_fasta)
		tot_repo_info_str_pdb.update(repo_info_str_pdb)
	fasta_repos = (tot_repo_info_seq_fasta, tot_repo_info_str_fasta)

	table = make_new_table(locations, fasta_repos, external_filename)

	return table
