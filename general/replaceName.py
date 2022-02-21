#!/usr/bin/python

import sys
import os

def replaceName(old, new, files):
    if '.' not in old:
        print("Phrase to replace does not appear to be a filename. Please include the file extension")
        sys.exit()
    for infile in files:
        if (infile[-2:] == '.o') or (infile[-2:] == '.e') or ('config' in infile):
            continue
        else: 
            print("Processing: ", infile)
            with open(infile, 'r') as file :
                    filedata = file.read()
            print('Occurences replaced: ', filedata.count(old))
            filedata = filedata.replace(old, new)

            #rewrite file with new data
            with open(infile, 'w') as file:
                file.write(filedata)

            if infile == (old):
                os.rename(infile, (new))

replaceName(sys.argv[1], sys.argv[2], sys.argv[3])