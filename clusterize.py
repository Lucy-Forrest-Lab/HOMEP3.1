# Name: clusterize.py
# Language: python3
# Libraries: 
# Description: Clusterizes the structures of HOMEP in structural families
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


def clusterize(locations, database, table, totaln_filename, HOMEP_filename, seqid_thr, tmscore_thr):
	obj_listofsets = []
	fam_listofsets = []


	if not table:
		table = {}
		totaln_file = open(locations['FSYS']['mainpath'] + totaln_filename, 'r')
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
			table[fields[0]][fields[1]][fields[2]][fields[3]] = (float(fields[4]), float(fields[5]), float(fields[6]), fields[7])

	structsets = []
	for ss in table:
		for superfamily in sorted(list(table[ss].keys)):
			structset = set()
			for struct in table[ss][superfamily]:
				structsets.append(frozenset(list(table[ss][superfamily][struct].keys())))

	for structset_1 in structsets:
		found = False
		for structset_2 in structsets:
			if structset_1 != structset_2 and structset_1 & structset_2:
				raise NameError("ERROR: intersection between two Superfamilies", structset_1 & structset_2)
			elif structset_1 == structset_2 and found:
				raise NameError("ERROR: more than one copy of the same structset")
			if structset_1 == structset_2:
				found = True

	HOMEP_file = open(locations['FSYS']['mainpath'] + HOMEP_filename, 'w')
	HOMEP_library = {}
	for ss in table:
		for sf in sorted(list(table[ss].keys)):
			superfamily_name = ss + str(int(sf))	
	
			obj_listofsets = []
			fam_listofsets = []
			for s1 in table[ss][sf]:
				for s2 in table[ss][sf][s1]:
					# Warning: there are cases in which seqid > seqid_thr but tmscore < tmscore_thr.
					# Usually they happen when one chain is a subsequence of the other one (in which case seqid ~ 1).
					# Nonetheless, we do not want to consider these cases as into the same object a priori.
					if table[ss][sf][s1][s2][1] > tmscore_thr:
						fam_listofsets.append(frozenset([s1, s2]))
					if table[ss][sf][s1][s2][0] > seqid_thr:
						obj_listofsets.append(frozenset([s1, s2]))

			obj_sets = merge(obj_listofsets)
			fam_sets = merge(fam_listofsets)

#		print("Superfamily #{0}".format(superfamily_name))
#		for obj in obj_sets:
#			print("Object")
#			for struct in obj:
#				print(struct)
	
			fam_obj_sets = []
			for nfam in range(len(fam_sets)):
				fam_obj_sets.append(set([]))
				for struct in fam_sets[nfam]:
					in_obj = False
					for nobj in range(len(obj_sets)):
						if struct in obj_sets[nobj]:
							fam_obj_sets[-1].add(obj_sets[nobj])
							in_obj = True
							break
					if not in_obj:
						fam_obj_sets[-1].add(frozenset([struct]))
	
			fl = []
			ol = []
			cf = 0
			for fam in fam_obj_sets:
				fl.append((fam, len(fam), cf))
				for obj in fam:
					ol.append((obj, len(obj), cf))
				cf += 1
	
			cf = 0
			co = 0
			HOMEP_file.write("Superfamily #{0}\n".format(superfamily_name))
			for fam in sorted(fl, key=lambda x : -x[1]):
				HOMEP_file.write("\tFamily #{0}\n".format(cf))
				for obj in sorted(ol, key=lambda y : -y[1]):
					if fam[2] != obj[2]:
						continue
					HOMEP_file.write("\t\tObject #{0}\n".format(co))
					for struct in obj[0]:
						HOMEP_file.write("\t\t\t{0}\t{1}\n".format(struct, database[struct][1]['FROM_PDB']['TITLE']))
					co += 1
				cf += 1
			HOMEP_library[superfamily_name] = fam_obj_sets
	HOMEP_file.close()
	return HOMEP_library
