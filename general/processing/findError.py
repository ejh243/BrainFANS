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
            outfile.write((re.findall(r"_(\d*)", file))[0]+'\t')
            with open(file+'.o', 'r') as output:
                for line in output:
                    sample = line.split(':  ')
                    if 'Current sample' in sample[0]:
                        outfile.write((sample[1].rstrip().lstrip())+'\t')
                    elif 'EXIT' in sample[0]:
                        outfile.write(sample[1].rstrip().lstrip()+'\t')
            with open(file+'.o', 'r') as file:
                content=file.read()
                if (error in content):
                        outfile.write(error+'\n')
                else:
                    outfile.write('\n')
                    
find_err(sys.argv[1], sys.argv[2:])