sampleName=$1
f1=$2
f2=$3

cd ${ALIGNEDDIR}
echo "Converting to genomedata format:" $sampleName

bamfile=$( find . -name ${sampleName}'*.bam' )

bedtools bamtobed -i ${bamfile} | genomedata-load -t stdin -s ${f1} -s ${f2} ${sampleName}.gnmdata