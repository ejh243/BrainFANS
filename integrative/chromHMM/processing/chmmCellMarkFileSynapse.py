#!/usr/bin/python

import sys

def CMBF_Synapse(file, mark):
	'''Function which takes the metadata provided by Synapse and outputs the cellMarkFileTable.txt with 
	case and control files required for the binarisation step of ChromHMM. 

	INPUT: 	synapse metadata file in comma-delimited format including the fastq file names
			the methylation mark of your files of interest
	OUTPUT: a txt file in the following format (tsv), named cellMarkFileTable.txt
		cell1 	mark1 	case1.bam 	control1.bam
		cell2	mark2	case2.bam 	control2.bam
	'''
	dic = {'oligodendrocyte': 'N-', 'GLUtamatergic neurons': 'N+', 'GABAergic neurons': 'N+'}
	case = []
	control=[]
	with open(file, 'r') as file:
		outfile = open('cellMarkFileTable.txt', 'wt')
		head = file.readline()
		for line in file:
			line = line.rstrip()
			line = line.replace('"', '')
			line = line.split(',') #split by column
			if line[10]=='input':
				control.append(line[4][:-12]+'_depDup_q30.bam')
			elif line[10] in mark:
				newline = dic[line[14]]+'\t'+line[10]+'\t'+line[4][:-12]+'_depDup_q30.bam'+'\t'
				if newline not in case:
					case.append(newline)
		control = list(set(control)) #remove duplicates generated from reads 1 and 2 
		#sort so that the lists are equivalent (by subject number and cell type)
		for z in range(2,0,-1):
			case.sort(key = lambda x: x.split(".")[z])
			control.sort(key = lambda x: x.split(".")[z])
		#write to outfile
		for x in range(0, len(case)):
			outfile.write(case[x]+control[x]+'\n')
	outfile.close()

CMBF_Synapse(sys.argv[1], sys.argv[2])