sampleName=$1
f1=$2
f2=$3

cd ${ALIGNEDDIR}
echo "Converting to genomedata format:" $sampleName

bamfile=$( find . -name ${sampleName}'*.bam' )

if [ ! "$(ls ${sampleName}.bed)" ]
then
	bedtools bamtobed -i ${bamfile} > ${sampleName}.bed
fi

genomedata-load -t ${sampleName}.bed -s ${f1} -s ${f2} ${sampleName}.gnmdata

if [[ $? == 0 ]]
then
	rm ${sampleName}.bed
	echo 'object generated'
fi
