## Performs fraction/group level peak calling using MACS with 1. shifted tagAlign files 
## and 2. paired end BAM files
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./ATACSeq/preprocessing/fractionPeaks.sh 

# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2

## INPUT
# array of shifted tag align files


## OUTPUT
# ${PEAKDIR}/MACS/ShiftedTagAlign/*.broadPeak (and other macs output)
# ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt

cd ${ALIGNEDDIR}

SAMPLEFILE=$1

QC1SAMPLES=($(awk '{print $1}' ${SAMPLEFILE} ))
SUBSET=( "${QC1SAMPLES[@]/%/.tn5.tagAlign.gz}" )

echo 'Files for peak-calling are:' ${SUBSET[@]}

macs2 callpeak -t ${SUBSET[@]} --outdir ${PEAKGROUP} -n ${GROUP}_atac -f BED -g 2.9e9 -q 5e-2 --keep-dup all  --shift 100 --extsize 147 --nomodel --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKGROUP}/${GROUP}_atac_peaks.broadPeak -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
	| grep -P 'chr[\dXY]+[ \t]' > ${PEAKGROUP}/${GROUP}_atac.broadPeak.filt

echo 'Finished calling peaks in common'