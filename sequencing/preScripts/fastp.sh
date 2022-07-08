## Performs trimming of paired fastq files fastp for a single sample

## EXECUTION
# sh ./fastp.sh <fastq file>
# where 
# <fastq file> is the path to the "R1" fastq files which are expected to be compressed, and have either r1 or R1 in the filename
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# RAWDATADIR, TRIMDIR

## REQUIRES the following software
# fastp

## INPUT
# 2 fastq files

## OUTPUT
# trimmed file
# trimming report 

sampleName=$1

echo
echo "Starting trimming on" ${sampleName} "at: "
date -u

f1=$(basename $2)
f2=$(basename $3)

cd ${RAWDATADIR}

## create output filenames
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}
  
mkdir -p ${TRIMDIR}/fastp_reports/
     
echo "Looking for trimmed files in" ${TRIMDIR}

if [ ! -s ${TRIMDIR}/fastp_reports/${sampleName}_fastp.json ]    
then
   echo "Trimmed fastq not found so trimming"
   fastp --detect_adapter_for_pe --length_required=27 --thread=8 --in1=${f1} --in2=${f2} --out1=${TRIMDIR}/${outf1} --out2=${TRIMDIR}/${outf2} --html=${TRIMDIR}/fastp_reports/${sampleName}_fastp.html --json=${TRIMDIR}/fastp_reports/${sampleName}_fastp.json
fi

if [[ $? == 0 ]]
   then echo "fastp complete"
fi
