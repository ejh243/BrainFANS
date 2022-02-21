#!/usr/bin/python

import sys
import re


def find_err(error, files):
    '''Function that takes the list of sample ids and searches a log file list for sample id and user-defined error messsage.
    
    INPUT: 
        error message to search for (provided on the command line)
        list of files to search
    OUTPUT:
        .log tsv file with sample ids in the first column, exit code in the second and error message in the third
    '''
    rgx = re.findall(r"(\d+)_", files[0])
    #create outfile name from first file
    outfile = (rgx[0]+'.log')
    with open(outfile, 'wt') as outfile:
        for file in files:
            with open(file+'.o', 'r') as output:
                for line in output:
                    if line[0:8] == 'Current ':
                        outfile.write((line[17:]).rstrip()+'\t')
                    elif line[0:8] == 'EXITCODE':
                        outfile.write(line[10:]+'\t')
            with open(file+'.e', 'r') as file:
                content=file.read()
                if (error in content):
                        outfile.write(error+'\n')
                else:
                    outfile.write('\n')

find_err(sys.argv[1], sys.argv[2], sys.argv[3:])