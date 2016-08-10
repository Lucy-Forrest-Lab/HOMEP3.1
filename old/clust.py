import sys

ftable = open(sys.argv[1], 'r')
text = ftable.read().split("\n")

clust = []
for line in text:
#	print(line)
	fields = line.split()
	if len(fields) != 4:
		continue
	found = False
	for i in range(0, len(clust)):
		if fields[0] in clust[i]:
			if (not (fields[1] in clust[i]) and float(fields[3])>0.85):
				clust[i].append(fields[1])
			found = True
			break
	if not found:
		clust.append([])
		clust[len(clust)-1].append(fields[0])
		if float(fields[3])>0.85:
			clust[len(clust)].append(fields[1])
#	print(clust)

#print(clust)
#print(len(clust), sum([len(i) for i in clust]))

cfile = open(sys.argv[2], 'r')
text = cfile.read().split("\n")

cl2cl = {}
for line in text:
	fields = line.split()
	if len(fields) > 0:
		if fields[0] == "group:":
			ng = int(fields[1])
		else:
			for i in range(0, len(clust)):
				if fields[0] in clust[i]:
					cl2cl[i] = ng
#print(cl2cl) 

for i in range(1, ng+1):
	print("Family #{0}".format(i))
	no = 0
	for j in cl2cl.keys():
		if cl2cl[j] == i:
			no += 1
			print("\tObject #{0}".format(no))
			for k in clust[j]:
				print("\t\t{0}".format(k))
