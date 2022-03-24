## script to process a single sample
## Requires a raw (r1) fastq file provided on the command line 
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	
## excludes duplicates, mt reads, only keeps properly paired reads
## shifts read prior to peak calling

sampleName=$1

echo "Changing to trimmed directory" $TRIMDIR
cd ${TRIMDIR}

## need trimmed files for alignment

## filename format changed by trimgalore to be val_1 and val_2 (not _trimmed)
f1=$(basename $(ls ${TRIMDIR}/${sampleName}*_val_1*.f*)) 
f2=$(basename $(ls ${TRIMDIR}/${sampleName}*_val_2*.f*)) 


echo "Found trimmed files:"
echo ${f1}
echo ${f2}

## checks if last file exists and is not empty

echo "Running alignment for" ${sampleName}
date -u	
## alignment
bismark --genome ${REFGENOME} -o ${ALIGNEDDIR} -1 $f1 -2 $f2

if [[ $? == 0 ]]
	then echo "Alignment and post filtering complete"
	date -u
fi