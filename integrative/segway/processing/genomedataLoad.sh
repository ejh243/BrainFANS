sampleName=$1
f=$2

echo
echo 'Starting creating genomedata object of' $sampleName 'at: '
date -u 

echo "Output directory is:" $LOADDIR

#create genomedata object

echo "Load genomedata sequence"
genomedata-load-seq --verbose ${LOADDIR}/${sampleName}.gnmdat ${REFGENOME}  

echo "Open genomedata trackfile"
genomedata-open-data --verbose ${LOADDIR}/${sampleName}.gnmdat --tracknames ${sampleName}

echo "Load genomedata trackfile"
genomedata-load-data --verbose ${LOADDIR}/${sampleName}.gnmdat ${sampleName} < ${LOADDIR}/${f} 
genomedata-close-data --verbose ${LOADDIR}/${sampleName}.gnmdat


if [[ $? == 0 ]]
then
	echo 'Object generated'
fi
