## calculates sequencing qc metrics with fastqc
## performs trimming of fastq files for a single sample
## expects paired end data where the two fastq files have either r1, r2,R1 or R2
## assumes sample name is located in filename directly before the specification of read type
## should accept any fastq suffix (e.g. .fastq, .fq) but compressed
## designed to be used with SLURM batch array jobs but could be used to process a single sample
## requires a (r1) fastq file provided on the command line 


f=$1 
cd ${DATADIRPE}
## extract sample names
FOLDER=$(dirname ${f})
f1=$(basename $f)
sampleName=${f1%_[rR]*}
echo "Processing" ${sampleName}


## create filename for paired fastq file
f2=$(basename $(ls ${FOLDER}/${sampleName}*[rR]2*q.gz))

## name output files
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}
  
## run fastqc
echo "Running FASTQC on"
echo ${f1}
echo ${f2}

echo "Output written to " ${FASTQCDIR}

cd ${FOLDER}  
fastqc ${f1} ${f2} -o ${FASTQCDIR}

## run fastp if not already run
echo "Running FASTP"

mkdir -p ${FOLDERTRIM}/fastp_reports/
 
echo "Looking for trimmed files in" ${FOLDERTRIM}

if [ ! -s ${FOLDERTRIM}/fastp_reports/${sampleName}_fastp.json ]	
then
  echo "Trimmed fastq not found so trimming"
  ## trim adapters only do not trim based on quality
  fastp --detect_adapter_for_pe --length_required=27 --thread=8 --in1=${f1} --in2=${f2} --out1=${FOLDERTRIM}/${outf1} --out2=${FOLDERTRIM}/${outf2} --html=${FOLDERTRIM}/fastp_reports/${sampleName}_fastp.html --json=${FOLDERTRIM}/fastp_reports/${sampleName}_fastp.json
fi

echo "FASTQC and fastp complete"
