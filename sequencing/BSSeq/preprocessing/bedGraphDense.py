#!/usr/bin/python

import sys

def dense(infile, outfile):
	'''Function that takes a bedgraph file as input and makes the chr start and end 
	positions continuous if line1:col3 == line2:col2.
	e.g.:

	in.bg
	chr1    10093   10094
	chr1    10094   10095
	chr1    10095   10096
	chr1    10096   10097
	chr1    10097   10098 
	chr1	10106	10116	

	out.bg
	chr1	10093	10098
	chr1	10106	10116

	INPUT: 
	    bedgraph file 
	OUTPUT:
	    bedgraph file with dense chromosome start end positions.
	'''
	infile=open(infile, 'r')
	outfile=open(outfile, 'wt')
	starts=[]
	lines=[infile.readline().rstrip()]
	for line in infile:
	    lines.append(line.rstrip())
	    line1=lines[0].split('\t')
	    line2=lines[1].split('\t')
	    starts.append(line1[1]) #add the chr start position
	    if line1[2] != line2[1]:        
	        outfile.write(line1[0]+'\t'+starts[0]+'\t'+line1[2]+'\n')
	        starts=[]
	    lines.pop(0)

	infile.close()
	outfile.close()

dense(sys.argv[1], sys.argv[2])