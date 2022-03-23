sampleName=$1
f1=$2
f2=$3

cd ${ALIGNEDDIR}
echo "Converting to genomedata format:" $sampleName

bamfile=$(basename $( find . -name ${sampleName}'*.bam' ))

if [ ! "$(ls ${sampleName}.bed)" ] || [ ! "$(ls ${sampleName}.tagAlign.gz)" ]
then
	bedtools bamtobed -i ${bamfile} > ${sampleName}.bed
fi

genomedata-load -s ${f1} -s ${f2} -t ${sampleName}.tagAlign.gz ${sampleName}.gnmdata

if [[ $? == 0 ]]
then
	rm ${sampleName}.bed
	echo 'Object generated'
fi
