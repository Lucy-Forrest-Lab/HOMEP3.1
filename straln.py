# Name: straln.py
# Language: python3
# Libraries: 
# Description: Runs the structural alignments needed for HOMEP
# Author: Edoardo Sarti
# Date: Aug 17 2016

import os, sys, multiprocessing, subprocess, re, time

def FrTMjob(data):
	ntm, maindir, checkpoint = data
	if not os.path.exists(maindir + '/' + ntm + '/alignments'):
		os.mkdir(maindir + '/' + ntm + '/alignments')
	if not os.path.exists(maindir + '/' + ntm + '/alignments/fasta'):
		os.mkdir(maindir + '/' + ntm + '/alignments/fasta')
	if not os.path.exists(maindir + '/' + ntm + '/alignments/str_alns'):
		os.mkdir(maindir + '/' + ntm + '/alignments/str_alns')
	if not os.path.exists(maindir + '/' + ntm + '/structures'):
		raise NameError("ERROR: The folder {0} is badly formatted and does not contain a structures/ subfolder.\n".format(maindir + '/' + ntm + '/') +
		                "       Please create one and fill it with all and only the appropriate pdb chains.")
	if os.path.exists(maindir + '/' + ntm + '/struct_codes.dat'):
		structcodesfile = open(maindir + '/' + ntm + '/struct_codes.dat', 'r')
		text = structcodesfile.read().split('\n')
		name2code = {}
		code2name = {}
		for line in text:
			fields = line.split()
			if len(fields) == 0:
				continue
			name2code[fields[1]] = [x for x in fields[0].split('.')]
			code2name[fields[0].split('.')[1]] = fields[1]
	else:
		raise NameError("ERROR: The folder {0} is badly formatted and does not contain a struct_codes.dat file.\n".format(maindir + '/' + ntm) +
		                "       Please generate it. It must contain all and only the names of the pdb chains in the structures/ subfolder"+
		                      " and each name must be associated with the correct structure code SC.\n" +
		                "       The format must be: <SC>\\t\\t<CHAIN NAME (XXXX_Y, no .pdb extension)>")
	for chain in name2code.keys():
		if not os.path.exists(maindir + '/' + ntm + '/structures/' + chain + '.pdb'):
			raise NameError("ERROR: The file {0} corresponding to Structure Code {1}".format(chain + '.pdb', name2code[chain]) +
			                      " was not found in the structures/ subfolder.")
	for struct in os.listdir(maindir + '/' + ntm + '/structures/'):
		if not struct[:6] in name2code:
			raise NameError("ERROR: The file {0} found in the structures/".format(struct) +
			                      " subfolder is not present in the struct_code.dat file.")
	if len(os.listdir(maindir + '/' + ntm + '/structures/')) < 2:
		return

	for chain_1 in [code2name[x] for x in sorted(code2name.keys())]:
#		file_1 = maindir + '/' + ntm + '/structures/' + chain_1 + '.pdb'
		file_1 = chain_1 + '.pdb'
		if os.path.exists(maindir + '/' + ntm + '/alignments/fasta/seq_'+chain_1+'.dat') and os.path.exists(maindir + '/' + ntm + '/alignments/str_alns/str_'+chain_1+'.dat'):
			continue
		if not os.path.exists(maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/'):
			os.mkdir(maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/')
		if not os.path.exists(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/'):
			os.mkdir(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/')
		for chain_2 in [code2name[x] for x in sorted(code2name.keys())]:
			if chain_1 == chain_2:
				continue
			print("#td "+ntm+"\t\tchain_1 "+chain_1+"\t\tchain_2 "+chain_2)
#			file_2 = maindir + '/' + ntm + '/structures/' + chain_2 + '.pdb'
			file_2 = chain_2 + '.pdb'
#			FTA_str_output = maindir + '/' + ntm + '/alignments/str_alns/tmp_' + chain_1 + '/' + chain_1 + '_' + chain_2 + '.tmp'
			FTA_str_output = name2code[chain_1][1] + '_' + name2code[chain_2][1] + '.tmp'
			FTA_seq_output = maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/' + name2code[chain_1][1] + '_' + name2code[chain_2][1] + '.tmp'

			FTA_stdout_file = open(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/aln_' + name2code[chain_1][1] + '_' + name2code[chain_2][1] + '.tmp', 'w')
			fnull = open(os.devnull, 'w')
			p = subprocess.Popen(['/v/apps/csb/frtmalign/frtmalign.exe', file_1, file_2, '-o', FTA_str_output], stdout=FTA_stdout_file, stderr=fnull, cwd=maindir+'/'+ntm+'/structures/')
			fnull.close()
			p.wait()
			FTA_stdout_file.close()
			os.rename(maindir + '/' + ntm + '/structures/' + FTA_str_output, maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/' + FTA_str_output)

			FTA_stdout_file = open(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/aln_' + name2code[chain_1][1] + '_' + name2code[chain_2][1] + '.tmp', 'r')
			text = FTA_stdout_file.read().split('\n')
			FTA_stdout_file.close()
			os.remove(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/aln_' + name2code[chain_1][1] + '_' + name2code[chain_2][1] + '.tmp')
			chkaln = -1000
			for nl in range(len(text)):
				if "Aligned length" in text[nl]:
					fields = re.split('=|,|\s',text[nl])
					fields = list(filter(None, fields))
#					print(fields)
					RMSD = float(fields[4])
					TMscore = float(fields[6])
				elif chkaln+1 == nl:
					seq_1 = text[nl]
				elif chkaln+3 == nl:
					seq_2 = text[nl]
				elif "denotes the residue pairs of distance" in text[nl]:
					chkaln = nl
			tmpseq_file = open(FTA_seq_output, 'w')
			tmpseq_file.write(">" + chain_1 + "\n" + seq_1.replace('\x00', '') + "\n>" + chain_2 + "\n" + seq_2.replace('\x00', '') + "\n\nRMSD\t{0:.2f}\nTM-score\t{1:.5f}\n\n".format(RMSD, TMscore))
			tmpseq_file.close()

		str_file = open(maindir + '/' + ntm + '/alignments/str_alns/str_' + chain_1 + '.dat', 'w')
		for tmp_filename in sorted(os.listdir(maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/')):
			chain_2_code = re.split('_|\.', tmp_filename)[-2]
#			print(chain_1, "chain_2_code "+chain_2_code, "tmp_filename "+tmp_filename, name2code)
			str_file.write("BEGIN \nCHAIN_1: " + chain_1  + "\nCHAIN_2: " +  code2name[chain_2_code] +
			               "\nSequence Alignment Code (SAC): " + name2code[chain_1][0] + 
			               "." + name2code[chain_1][1] + "." + chain_2_code + "\n")
			tmp_file = open(maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/' + tmp_filename)
			text = tmp_file.read().split('\n')
			for line in text:
				str_file.write(line+'\n')
			str_file.write("END\n\n\n")
			os.remove(maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/' + tmp_filename)
			tmp_file.close()
		time.sleep(1)
		os.rmdir(maindir + '/' + ntm + '/alignments/str_alns/tmp_' + name2code[chain_1][1] + '/')
		str_file.close()

		seq_file = open(maindir + '/' + ntm + '/alignments/fasta/seq_' + chain_1 + '.dat', 'w')
		for tmp_filename in sorted(os.listdir(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/')):
			chain_2_code = re.split('_|\.', tmp_filename)[-2]
#			print(chain_1, "chain_2_code "+chain_2_code, "tmp_filename "+tmp_filename, name2code)
			seq_file.write("BEGIN \nCHAIN_1: " + chain_1  + "\nCHAIN_2: " +  code2name[chain_2_code] +
			               "\nSequence Alignment Code (SAC): " + name2code[chain_1][0] + 
			               "." + name2code[chain_1][1] + "." + chain_2_code + "\n")
			FTA_seq_output = maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/' + name2code[chain_1][1] + '_' + chain_2_code + '.tmp'
			tmp_file = open(FTA_seq_output, 'r')
			text = tmp_file.read().split('\n')
			for line in text:
				seq_file.write(line+'\n')
			seq_file.write("END\n\n\n")
			os.remove(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/' + tmp_filename)
			tmp_file.close()
		time.sleep(5)
#		print(os.listdir(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/'))
		os.rmdir(maindir + '/' + ntm + '/alignments/fasta/tmp_' + name2code[chain_1][1] + '/')
		seq_file.close()


def structure_alignment(locations, aligner_path, np):
	already_processed = []
	hidden_repository_filename = locations['mainpath'] + '.structure_alignments.dat'
	if os.path.exists(hidden_repository_filename):
		hidden_repository_file = open(hidden_repository_filename, 'r')
		text = hidden_repository_file.read().split("\n")
		for line in text:
			fields = line.split()
			already_processed.append((fields[2], fields[3]))
	
	subdirs = []
	for i in os.listdir(str(sys.argv[1])):
		if re.match('^\d*$', str(i)):
			subdirs.append(int(i))
	superfamilies = [(str(i), maindir+'/', 0) for i in sorted(subdirs)]
	for sf in superfamilies:
		
		
	pass	
"""
if len(sys.argv) < 2:
        raise NameError("Usage: start_FrTM.py <filesystem_main_dir> [{<subdir_names>}]")
maindir = sys.argv[1]
if not os.path.exists(maindir):
	raise NameError("ERROR: Directory {0} does not exists.".format(maindir))

nsubdirs = len(sys.argv) - 2
if nsubdirs > 0:
	subdirs = []
	for i in range(0, nsubdirs):
		subdirs.append(int(sys.argv[2+i]))
else:
	subdirs = []
	for i in os.listdir(str(sys.argv[1])):
		if re.match('^\d*$', str(i)) and os.path.exists(maindir + '/' + str(i) + '/struct_codes.dat'):
			subdirs.append(int(i))

superfamilies = [(str(i), maindir+'/', 0) for i in sorted(subdirs)]
#print(superfamilies)

#for sf in superfamilies:
#	FrTMjob(sf)


#exit(1)

pool = multiprocessing.Pool(processes=4)
pool_outputs = pool.map(FrTMjob, superfamilies)
"""
