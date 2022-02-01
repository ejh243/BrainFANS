## Calculates sequencing qc metrics with fastqc for paired fastq files for a single sample

## EXECUTION
# sh ./fastqc.sh <fastq file>
# where 
# <fastq file> is the path to the "R1" fastq files which are expected to be compressed, and have either r1 or R1 in the filename
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

f=$1

## extract sample name from filename
f1=$(basename $f)
sampleName=${f1%_[rR][12]*} ## sample name is everything before either r1 or R1
echo "Processing" ${sampleName}

## create filename for paired fastq file
f2=$(basename $(ls ${RAWDATADIR}/${sampleName}*[rR]2*q.gz))

  ## run fastqc
echo "Running FASTQC on"
echo ${f1}
echo ${f2}

echo "Output written to " ${FASTQCDIR}

cd ${RAWDATADIR}  
fastqc ${f1} ${f2} -t 8 -o ${FASTQCDIR}

#If fastqc is finished on all files then run MultiQC 
if [ `expr ${SLURM_ARRAY_TASK_ID} + 1` = "${#FQFILES[@]}" ]
then 
  echo "Running MultiQC on all files"
  multiqc ./*_fastqc.zip -o ${RAWDATADIR}
fi

echo "FASTQC and MultiQC complete"
