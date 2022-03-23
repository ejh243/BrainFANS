## Calculates sequencing qc metrics with fastqc for paired fastq files for a single sample

## EXECUTION
# sh ./fastqc.sh <sample name> <fastq R1 file> <fastq R2 file>
# where 
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# RAWDATADIR, FASTQCDIR

## REQUIRES the following software
# fastqc, multiqc,

## INPUT
# 2 fastq files

## OUTPUT
# fastqc reports for each fastq file
# multiqc on total fastqc output


sampleName=$1
f1=$(basename $2)
f2=$(basename $3)

cd ${RAWDATADIR} 
 
## extract sample name from filename
echo "Processing" ${sampleName}

## run fastqc

echo "Running FASTQC on"
echo ${f1}
echo ${f2}

echo "Output written to " ${FASTQCDIR}

fastqc ${f1} ${f2} -t 8 -o ${FASTQCDIR}

if [[ $? == 0 ]]
	then echo "FASTQC complete"
fi
