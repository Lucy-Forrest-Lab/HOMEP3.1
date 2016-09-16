# Name: cluster.py
# Language: python3
# Libraries: 
# Description: Clusterizes the structures of HOMEP in structural foldilies
# Author: Edoardo Sarti
# Date: Aug 17 2016

def merge(lsts):
	sets = [set(lst) for lst in lsts if lst]
	merged = 1
	while merged:
		merged = 0
		results = []
		while sets:
			common, rest = sets[0], sets[1:]
			sets = []
			for x in rest:
				if x.isdisjoint(common):
					sets.append(x)
				else:
					merged = 1
					common |= x
			results.append(frozenset([i for i in common]))
		sets = results
	return sets


def cluster(options, locations, database, table):
	totaln_filename = options['output_tab']
	HOMEP_filename = options['output_homep']
	seqid_thr = float(options['subunit_thr'])
	tmscore_thr = float(options['cluster_thr'])

	sub_listofsets = []
	fold_listofsets = []

	if not table:
		table = {}
		totaln_file = open(locations['FSYSPATH']['main'] + totaln_filename, 'r')
		text = totaln_file.read().split("\n")
		for line in text:
			if not line:
				continue
			fields = line.split()
			if not fields[0] in table:
				table[fields[0]] = {}
			if not fields[1] in table[fields[0]]:
				table[fields[0]][fields[1]] = {}
			if not fields[2] in table[fields[0]][fields[1]]:
				table[fields[0]][fields[1]][fields[2]] = {}
			table[fields[0]][fields[1]][fields[2]][fields[3]] = (float(fields[4]), float(fields[5]), float(fields[6]), float(fields[7]))


	HOMEP_file = open(locations['FSYSPATH']['main'] + HOMEP_filename, 'w')
	HOMEP_library = {}
	for toptype in sorted(list(table.keys())):
		for top in sorted(list(table[toptype].keys()), key = lambda x: int(x)):
			topology_name = toptype + str(int(top))	

			sub_listofsets = []
			fold_listofsets = []
			for s1 in table[ss][sf]:
				for s2 in table[ss][sf][s1]:
					# Warning: there are cases in which seqid > seqid_thr but tmscore < tmscore_thr.
					# Usually they happen when one chain is a subsequence of the other one (in which case seqid ~ 1).
					# Nonetheless, we do not want to consider these cases as into the same subunit a priori.
					if table[ss][sf][s1][s2][2] > tmscore_thr:
						fold_listofsets.append(frozenset([s1, s2]))
					if table[ss][sf][s1][s2][0] > seqid_thr:
						sub_listofsets.append(frozenset([s1, s2]))			

			# Creates closed component stes
			sub_sets = merge(sub_listofsets)
			fold_sets = merge(fold_listofsets)

			fold_sub_sets = []
			for nfold in range(len(fold_sets)):
				fold_sub_sets.append(set([]))
				# Search through all fold clusters
				for struct in fold_sets[nfold]:
					# The structure is not a priori in any seq cluster
					in_sub = False
					# Search through all sequence clusters
					for nsub in range(len(sub_sets)):
						# If a structure is in that seq cluster,
						# add that entire seq cluster to the last fold set
						if struct in sub_sets[nsub]:
							fold_sub_sets[-1].add(sub_sets[nsub])
							in_sub = True
							break
					# If struct was not in any seq cluster, add only that struct
					# to the last fold set
					if not in_sub:
						fold_sub_sets[-1].add(frozenset([struct]))

			# Create list of folds and list of subunits, each fold and subunit reports the number of elements it contains
			fl = []
			sl = []
			cf = 0  # Index to relate a fold to the subunits it contains
			for fold in fold_sub_sets:
				fl.append((fold, len(fold), cf))
				for sub in fold:
					sl.append((sub, len(sub), cf))
				cf += 1

			# Write the HOMEP file
			cf = 0
			cs = 0
			tmlist = []  # List to keep track of all classified stuctures
			HOMEP_file.write("Topology #{0}\n".format(topology_name))
			for fold in sorted(fl, key=lambda x : -x[1]):
				HOMEP_file.write("\tFold #{0}\n".format(cf))
				for sub in sorted(sl, key=lambda y : -y[1]):
					# If that subunit does not belong to that fold, jump
					if fold[2] != sub[2]:
						continue
					HOMEP_file.write("\t\tSubunit #{0}\n".format(cs))
					for struct in sub[0]:
						title = ''
						if 'TITLE' in database[struct[:4]][1]['FROM_PDB']:
							title = database[struct[:4]][1]['FROM_PDB']['TITLE']
						HOMEP_file.write("\t\t\t{0}\t{1}\n".format(struct, title))
						tmlist.append(struct)
					cs += 1
				cf += 1

			print(tmlist)

			# If there are other unclassified structures (not related to any other structure by any common fold or subunit),
			# for each of them a separate fold is created (super-singleton)
			for s in list(database.keys()):
				for c in list(database[s][1]['CHAIN'].keys()):
					struct = s+'_'+c
					if int(database[s][1]['CHAIN'][c][0]['NUM_TM']) == int(sf) and database[s][1]['CHAIN'][c][0]['TYPE'] == ss and struct not in tmlist:
						print(struct)
						HOMEP_file.write("\tFold #{0}\n".format(cf))
						HOMEP_file.write("\t\tSubunit #{0}\n".format(cs))
						title = ''
						if 'TITLE' in database[s][1]['FROM_PDB']:
							title = database[struct[:4]][1]['FROM_PDB']['TITLE']
						HOMEP_file.write("\t\t\t{0}\t{1}\n".format(struct, title))
						cs += 1
						cf += 1

			HOMEP_library[topology_name] = fold_sub_sets
	HOMEP_file.close()
	return HOMEP_library
