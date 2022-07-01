## script to process a single sample
## Requires sample name provided on the command line 
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	


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

echo "Running alignment for" ${sampleName}
date -u	

# alignment
bismark --genome ${REFGENOME} -o ${ALIGNEDDIR} -1 $f1 -2 $f2 --basename ${sampleName} --parallel

# align to spike-in genome
bismark --genome ${REFSPIKE} -o ${ALIGNEDDIR}/spikeAlignments -1 $f1 -2 $f2 --basename ${sampleName}.spike --parallel

# deduplicate for WGBS libraries
cd ${ALIGNEDDIR}
echo 'Deduplicating'
deduplicate_bismark --bam ${sampleName}*pe.bam -p

# extract context-dependent methylation
cd ${ALIGNEDDIR}
echo 'Extracting methylation'
bismark_methylation_extractor -p ${sampleName}*deduplicated.bam -o ${METHYLDIR} 

if [[ $? == 0 ]]
	then echo "Alignment and methylation extraction complete"
	date -u
fi