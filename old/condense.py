# Name: condense.py
# Language: python3
# Libraries: sys, os, re
# Description: Organizes alignment files in a compact file system
# Author: Edoardo Sarti
# Date: Aug 04 2016

import sys, os, re

if len(sys.argv) < 3:
	raise NameError("Usage: condense.py <new_filesystem_main_dir> <old_filesystem_main_dir> [{<subdir_names>}]")
if sys.argv[1] == sys.argv[2]:
	raise NameError("ERROR: New and old filesystems cannot have the same main directory")

nsubdirs = len(sys.argv) - 3
if nsubdirs > 0:
	subdirs = []
	for i in range(0, nsubdirs):
		subdirs.append(int(sys.argv[3+i]))
else:
	subdirs = []
	for i in os.listdir(str(sys.argv[2])):
		if re.match('^\d*$', str(i)):
			subdirs.append(int(i))

subdirs = [str(i) for i in sorted(subdirs)]

print(subdirs)

for i in subdirs:
	if not os.path.exists(sys.argv[1]+'/'+i):
		os.mkdir(sys.argv[1]+'/'+i)
	if not os.path.exists(sys.argv[1]+'/'+i+'/alignments'):
		os.mkdir(sys.argv[1]+'/'+i+'/alignments')
	if not os.path.exists(sys.argv[1]+'/'+i+'/alignments/fasta'):
		os.mkdir(sys.argv[1]+'/'+i+'/alignments/fasta')
	if not os.path.exists(sys.argv[1]+'/'+i+'/alignments/str_alns'):
		os.mkdir(sys.argv[1]+'/'+i+'/alignments/str_alns')
	structures = []
	for j in os.listdir(sys.argv[2]+'/'+i):
		if re.match('^\S*pdb\s*$', str(j)):
			structures.append(str(j)[:6])

	counter = 0
	count_j1 = 0
	codefile = open(sys.argv[1]+'/'+i+'/struct_codes.dat', 'w')
	for j1 in structures:
		count_j1 += 1
		count_j2 = 0
		new_seqfile = open(sys.argv[1]+'/'+i+'/alignments/fasta/seq_'+j1+'.dat', 'w')
		new_strfile = open(sys.argv[1]+'/'+i+'/alignments/str_alns/str_'+j1+'.dat', 'w')
		codefile.write('a'+i.zfill(3)+'.'+str(count_j1).zfill(6)+"    "+j1+"\n")
		for j2 in structures:
			count_j2 += 1
			if j1 == j2:
				continue
			counter += 1
			SAC = 'a'+i.zfill(3)+'.'+str(count_j1).zfill(6)+'.'+str(count_j2).zfill(6)
#			print("{0}".format(i+'/alignments/fasta/'+j1+'_'+j2+'.fasta.txt'))
			seqfile = open(sys.argv[2]+'/'+i+'/alignments/fasta/'+j1+'_'+j2+'.fasta.txt', 'r')
			strfile = open(sys.argv[2]+'/'+i+'/alignments/str_alns/'+j1+'_'+j2+'.sup', 'r')
			seqtext = seqfile.read().split('\n')
			strtext = strfile.read().split('\n')
			new_seqfile.write("BEGIN    "+str(counter)+"\n")
			new_seqfile.write("CHAIN_1: "+j1+"\nCHAIN_2: "+j2+"\nSequence Alignment Code (SAC): "+SAC+"\n")
			for l in seqtext:
				new_seqfile.write(l+"\n")
			new_seqfile.write("END\n\n\n")
			seqfile.close()
			new_strfile.write("BEGIN    "+str(counter)+"\n")
			new_strfile.write("CHAIN_1: "+j1+"\nCHAIN_2: "+j2+"\nSequence Alignment Code (SAC): "+SAC+"\n")
			for l in strtext:
				new_strfile.write(l+"\n")
			new_strfile.write("END\n\n\n")
			strfile.close()
		new_seqfile.close()
		new_strfile.close()
