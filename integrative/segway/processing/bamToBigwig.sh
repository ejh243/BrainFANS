sampleName=$1
f=$2

echo "Generating coverage file"

bamCoverage -b ${f} -o ${ALIGNEDDIR}/coverage/${sampleName}.bw -of 'bigwig' --region chr1 

if [[ $? == 0 ]]
then
	echo 'File generated'
fi