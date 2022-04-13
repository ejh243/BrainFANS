sampleName=$1
f=$2


cd ${INDATADIR}
echo "Converting to genomedata format:" $sampleName

#create genomedata object

echo 'Load genomedata sequence'
genomedata-load-seq ${sampleName}.gnmdat ${REFGENOME} --verbose 

echo 'Open genomedata trackfile'
genomedata-open-data ${sampleName}.gnmdat --tracknames high --verbose

echo 'Load genomedata trackfile'
genomedata-load-data ${sampleName}.gnmdat high < ${f} --verbose
genomedata-close-data --verbose ${sampleName}.gnmdat

if [[ $? == 0 ]]
then
	echo 'Object generated'
fi
