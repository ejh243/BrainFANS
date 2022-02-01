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


f=$1

cd ${RAWDATADIR}

## extract sample name from filename
f1=$(basename $f)
sampleName=${f1%_[rR][12]*} ## sample name is everything before either r1 or R1
echo "Processing" ${sampleName}

## create filename for paired fastq file
f2=$(basename $(ls ${RAWDATADIR}/${sampleName}*[rR]2*q.gz))

## create output filenames
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}
  
echo "Running FASTP"

mkdir -p ${TRIMDIR}/fastp_reports/
     
echo "Looking for trimmed files in" ${TRIMDIR}

if [ ! -s ${TRIMDIR}/fastp_reports/${sampleName}_fastp.json ]    
then
   echo "Trimmed fastq not found so trimming"
   fastp --detect_adapter_for_pe --length_required=27 --thread=8 --in1=${f1} --in2=${f2} --out1=${TRIMDIR}/${outf1} --out2=${TRIMDIR}/${outf2} --html=${TRIMDIR}/fastp_reports/${sampleName}_fastp.html --json=${TRIMDIR}/fastp_reports/${sampleName}_fastp.json
fi

echo "fastp complete"
