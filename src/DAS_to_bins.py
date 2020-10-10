#!/usr/bin/env python
### This script is used for extract contigs from a large contig file according to the DAS tools output (contig identifier VS bin name)
### jianshu.zhao@gatech.edu
## usage: python DAS_to_bin.py contigs.fa DASTool_scaffolds2bin.txt GWMC_1222_(output prefix)

import sys
file1=str(sys.argv[1])
file2=str(sys.argv[2])
prefix=str(sys.argv[3])
## read scaffolds2bin.txt
rows=""
dict2 = {}
with open(file2, "r") as f:
	for line in f:
		line = line.strip()
		rows = line.split("\t")
		dict2[">" + rows[0]] = rows[1]
### read contig file
dict1={}
table2 = []
rows1 = ""
with open(file1, "r") as f:
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			rows = line.split(" ")
			rows1 = line.strip()
			table2.append(rows1)
		else:
			table2.append(line)
for i in range(0,len(table2),2):
	dict1[table2[i]] =  table2[i+1]

for i in dict1:
	if str(i)[:str(i).index(" ")] in dict2:
		file_fh = prefix + str(dict2[str(i)[:str(i).index(" ")]]) + ".fasta"
		with open(file_fh, "a") as fh:
			fh.write(str(i) + "\n" + str(dict1[i]) + "\n")
