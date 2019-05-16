## using output of king and plink --miss identify worse performing duplicate for exclusion 
## this script is run as follows
## python ExcludeDuplicates.py <king output> <plink --miss output> <output file>


import sys


print "Reading in sample missingness from", sys.argv[2]

sampleMissing = file(sys.argv[2], "r")
sampleMissing = sampleMissing.readlines()

sMissing = {}
for line in sampleMissing:
    line = line.strip().split()
    sMissing[line[0]+"^"+line[1]] = line[4]
print "Missingness info read for ", str(len(sMissing)), "samples"	
print "\n\n"


print "Identifying duplicates from", sys.argv[1]

dupSamplesFile = file(sys.argv[1], "r")
dupSamplesFile = dupSamplesFile.readlines()

dupSamples = {} 
for line in dupSamplesFile:
	line = line.strip().split()
	sample1 = line[0] + "^" + line[1]
	sample2 = line[2] + "^" + line[3]
	
	## need to save both ways around in case 3 way+ duplicates
	if sample1 in dupSamples:
	   oldEntry = dupSamples[sample1]
	   newEntry = oldEntry + [sample2]
	   dupSamples[sample1] = newEntry
	else: 
	   dupSamples[sample1] = [sample2]
	   
	if sample2 in dupSamples:
	   oldEntry = dupSamples[sample2]
	   newEntry = oldEntry + [sample1]
	   dupSamples[sample2] = newEntry
	else: 
	   dupSamples[sample2] = [sample1]

## create unique list of duplicate samples
fullList = []
for item in dupSamples:
	allDups = dupSamples[item] + [item]
	allDups.sort()
	fullList = fullList + [";".join(allDups)]

uniqList = []
unique = [uniqList.append(x) for x in fullList if x not in uniqList]
print str(len(uniqList)),"groups of duplicate samples"
uniqList.sort()


print "Writing list of samples to exclude to", sys.argv[3]
output = file(sys.argv[3], "w")
## find sample with least missingness and exclude all others
for item in uniqList:
	samples = item.split(";")
	list_values = []
	
	for element in samples:
		if element in sMissing:
			list_values = list_values + [sMissing[element]]
			
		else:
			print "Couldn't find sample: ", element
	
	if len(list_values) != 0:
		indexToKeep = list_values.index(min(list_values))
		samples.remove(samples[indexToKeep])
		for each in samples:

			output.write("\t".join(each.split("^")) + "\n")

output.close()
	   

