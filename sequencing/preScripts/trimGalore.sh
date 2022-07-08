## performs trimming of fastq files for a single sample
## expects paired end data where the two fastq files are specified _1 and _2
## should accept any fastq suffix (e.g. .fastq, .fq) but compressed
## designed to be used with SLURM batch array jobs but could be used to process a single sample
## requires a (_1) fastq file provided on the command line 

sampleName=$1
echo
echo "Starting trimming on" ${sampleName} "at: "
date -u

f1=$(basename $2)
f2=$(basename $3)

## create output filenames
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}
  
## run trim galore if not already run
echo "Running Trim Galore"

cd ${RAWDATADIR}
mkdir -p ${TRIMDIR}/trimGaloreReports/

echo "Looking for trimmed files in" ${TRIMDIR}
if [ -z "$(ls -A ${TRIMDIR}/trimGaloreReports)" ]; 
then
  echo "Trimmed fastq not found so trimming"
  trim_galore --paired ${f1} ${f2} -o ${TRIMDIR} --clip_R2 2
  ##quality score default 20
  ##default min length 20
  ##clip_R2 used specifically for WGBS read 2 due to end repair reaction introducing methylation biases in PE
fi

if [[ $? == 0 ]]
then 
  mv ${sampleName}*report.txt trimGaloreReports
  echo 'Trim Galore complete'
fi