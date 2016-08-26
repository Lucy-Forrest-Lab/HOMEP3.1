# Name: opm.py
# Language: python3
# Libraries: 
# Description: Library of functions in order to run PPM and get/store OPM data
# Author: Antoniya Aleksandrova, Edoardo Sarti
# Date: Aug 15 2016

import os, shlex, multiprocessing

def TMdoms():
	main_dir = '/u/aa606/NIH_work/'
	list_file = open(main_dir+'ppm_results/alpha_list.txt', 'r')
	opm_data = {}
	for pdb in list_file.read().split('\n'):
		if not pdb:
			continue
		results = open(main_dir + 'ppm_results/data/' + pdb[0:4] + '_results', 'r')
		units = ''
		flag = 0
		for line in results:
			if 'possible transmembrane secondary structure segments' in line:
				fields=line.split()
				if fields[0] == '0':
					flag = 1
				break
		results.close()
		if flag == 0:
			datasub_file = open(main_dir+'ppm_results/data/' + pdb[0:4] + '_datasub1', 'r')
			pdbname = pdb[0:4].upper()
			opm_data[pdbname] = {}
			for line in datasub_file:
				if pdb[0:4] in line:
					seg = 0
					fields = line.split(';')
					chain = fields[1]
					seg = len(fields[3].split(','))
					opm_data[pdbname][chain] = seg
			datasub_file.close()
	list_file.close()
	return opm_data


def run_PPM(data):
	locations, struct = data
	main_dir = '/u/aa606/NIH_work/'
	ppm_path = ''

	ppminput_filename = ''
	ppminput_file = open(ppminput_filename, 'w')
	ppminput.write(' 0 in  '+struct_path+struct+'.pdb')
	p = subprocess.Popen(ppm_path, stdout=, stdin=ppminput_file)
	ppminput_file.close()


def run_CEsym(data):
	pass


def run_SYMD(data):


def extract_fields(field, sep):
	if type(sep) != str:
		raise NameError("ERROR: separator is not a string")
	if len(sep) > 1:
		raise NameError("ERROR: separator is not a single character")
	if type(field) != str:
		raise NameError("ERROR: only strings can be separated in strings")
	subfields = field.split(sep)
	if len(subfields) > 1:
		lst = []
		for nsf in range(len(subfields)):
			lst.append(subfields[nsf])
		return lst
	else:
		return field


def read_table(locations)
	table_filename = locations['OPT']['OPMTABLE']
	table_file = open(table_filename, 'r')
	text = table_file.read().split('\n')
	table_file.close()

	OPM_dict = {}
	attr = []
	line = text[0]
	fields = line.split()
	for f in fields[2:]:
		attr.append(f)
	for line in text[1:]:
		if not line:
			continue
		fields = shlex.split(line)
		subfields = fields[1].split(';')
		if len(subfields) == 1:
			struct = fields[0] + '_' + subfields[0]
			if struct not in OPM_dict:
				OPM_dict[struct] = {}
			for nf in range(len(fields[2:])):
				if ';' in fields[nf]:
					OPM_dict[attr[nf-2]] = extract_fields(fields[nf], ';')
				elif ',' in fields[nf]:
					OPM_dict[attr[nf-2]] = extract_fields(fields[nf], ',')
				else:
					OPM_dict[attr[nf-2]] = fields[nf]
		else:
			for nsf in range(subfields):
				struct = fields[0] + '_' + subfields[nsf]
				if struct not in OPM_dict:
					OPM_dict[struct] = {}
				for nf in range(len(fields[2:])):
					if ';' in fields[nf]:
						OPM_dict[attr[nf-2]] = extract_fields(fields[nf], ';')
					elif ',' in fields[nf]:
						OPM_dict[attr[nf-2]] = extract_fields(fields[nf], ',')
					else:
						OPM_dict[attr[nf-2]] = fields[nf]	
	return OPM_dict			


def OPM_retriever(locations, structlist, np):
	# If struct is not found in the table, run the scripts and update the table.
	# Retrieve information from the table.

	OPM_dict = read_table(locations)

	data = []
	for struct in structlist:
		if struct not in OPM_dict:
			data.append((locations, struct))

	pool = multiprocessing.Pool(processes=np)
	pool_outputs = pool.map(FrTMjob, data)

	return OPM_dict
