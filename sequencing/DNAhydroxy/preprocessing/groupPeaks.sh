## Performs fraction/group level peak calling using MACS 
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./DNAhydroxy/preprocessing/fractionPeaks.sh 

# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2

## INPUT
# array of bam files


## OUTPUT
# ${PEAKDIR}/MACS/subset

cd ${ALIGNEDDIR}

SAMPLEFILE=$1

QC1SAMPLES=($(awk '{print $1}' ${SAMPLEFILE} ))
SUBSET=( "${QC1SAMPLES[@]/%/*.deduplicated.bam}" )

CONTROL=($(grep "${GROUP}" ${SAMPLEFILE} | awk '{print $4}'))
CONTROL=( "${CONTROL[@]/%/*.deduplicated.bam}" )

mark=($(awk '{print $3}' $METADIR/samplesForGroupAnalysis.txt | sort -u))

echo 'Files are:' ${SUBSET[@]}
echo 'Controls are: ' ${CONTROL[@]}

macs2 callpeak -t ${SUBSET[@]} -c ${CONTROL[@]} --outdir ${PEAKGROUP} -f BAM -g 2.9e9 -n ${GROUP}_${mark} -q 1e-2 

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKGROUP}/${GROUP}_${mark}*peaks.narrowPeak -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
	| grep -P 'chr[\dXY]+[ \t]' > ${PEAKGROUP}/${GROUP}_${mark}.narrowPeak.filt