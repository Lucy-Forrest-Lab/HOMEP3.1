# Name: OPM_TMdoms.py
# Language: python3
# Libraries: 
# Description: Checks TM domains on OPM
# Author: Antoniya Aleksandrova, Edoardo Sarti
# Date: Aug 15 2016

import os

def OPM_TMdoms():
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
