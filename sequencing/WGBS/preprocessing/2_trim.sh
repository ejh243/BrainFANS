## performs trimming of fastq files for a single sample
## expects paired end data where the two fastq files are specified _1 and _2
## assumes sample name is located in filename directly before the specification of read type
## should accept any fastq suffix (e.g. .fastq, .fq) but compressed
## designed to be used with SLURM batch array jobs but could be used to process a single sample
## requires a (_1) fastq file provided on the command line 


f=$1 
cd ${DATADIR}
## extract sample names
FOLDER=$(dirname ${f})
f1=$(basename $f)
sampleName=${f1%_*1.*} ##rm [rR] add 1*
echo "Processing" ${sampleName}


## create filename for paired fastq file
f2=$(basename $(ls ${FOLDER}/${sampleName}*_*2.*q.gz)) ##rm [rR]

## name output files
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}

## run trim galore if not already run
echo "Running Trim Galore"

mkdir -p ${FOLDERTRIM}/trimGaloreReports/

cd $FOLDER 
echo "Looking for trimmed files in" ${FOLDERTRIM}
if [ -z "$(ls -A ${FOLDERTRIM}/trimGaloreReports)" ]; 
then
  echo "Trimmed fastq not found so trimming"
  trim_galore --paired ${f1} ${f2} -o ${FOLDERTRIM} --clip_R2 2
  ##quality score default 20
  ##default min length 20
  ##clip_R2 used specifically for WGBS read 2 due to end repair reaction introducing methylation biases in PE
fi

##run fastqc on the trimmed files

module purge
module load FastQC
module load MultiQC

cd $FOLDERTRIM
fastqc ${outf1} ${outf2} -t 8 -o ${FASTQCDIR}

#If all files are present then run MultiQC 
if [ `expr ${SLURM_ARRAY_TASK_ID} + 1` = "${#FQFILES[@]}" ]
then 
  echo "Running MultiQC on all files"
  multiqc ./*_fastqc.zip
fi