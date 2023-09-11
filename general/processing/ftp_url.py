#!/usr/bin/python

import sys

def ftp(infile):
    '''Function to use filereport from ENA to generate the ftp URLs for download
    INPUT: filereport from ENA
    OUTPUT: text list of ftp URLs'''
    with open(infile, 'r') as infile:
        with open('ftp_url_list.txt', 'wt') as outfile:
            header = infile.readline()
            for line in infile:
                line = line.split('\t')
                line = line[6].split(';')
                line = '\nftp://'.join(line)
                outfile.write('ftp://'+line+'\n')
				
ftp(sys.argv[1])