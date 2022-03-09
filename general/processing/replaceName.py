#!/usr/bin/python

import sys
import os

def replaceName(old, new, files):
    '''Function to replace all filenames within all scripts in a directory. 

    INPUT:  phrase to replace including the file extension
            phrase to replace with
            list of files to find an replace within 
    '''
    if not any(item in old for item in ['/', '.']):
        print("Phrase to replace does not appear to be a filename. Please include the file extension or directory / at the end")
        sys.exit()
    files=files[0].split('\n')
    for infile in files:
        if (infile[-2:] == '.o') or (infile[-2:] == '.e') or ('config' in infile):
            continue
        else: 
            print('Processing: '+infile)
            with open(infile, 'r') as file :
                    filedata = file.read()
            print('Occurrences replaced: '+str(filedata.count(old)))
            filedata = filedata.replace(old, new)

            #rewrite file with new data
            with open(infile, 'w') as file:
                file.write(filedata)

            if (os.path.basename(infile) == (old)):
                DIR=os.path.dirname(infile)
                print('Renaming: '+infile+' to '+(DIR+'/'+new))
                DIR=os.path.dirname(infile)
                os.rename(infile, (DIR+'/'+new))

            print('')

replaceName(sys.argv[1], sys.argv[2], sys.argv[3:])