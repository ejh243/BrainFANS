#!/bin/bash

## ================================================================================================##
##                           ATAC-seq pipeline STEP 5.2: Sex Chr Peak calling                      ||                  
## ================================================================================================##
## EXECUTION: sh ./subScripts/sexChrPeaks.sh                                                       ||
## - execute from scripts directory                                                                ||
##                                                                                                 ||
## DESCRIPTION: This script call peaks only on reads assigned to sex chromosomes merging all       ||
## samples, using MACS3 Single-end mode                                                            ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## ${PEAK_DIR}/ShiftedTagAlign/sexChr/chr*.broadPeak.filt                                          ||
## - output will be peaks called in sex chromosomes and not aligned to par regions                 ||
## ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/chr*.peakcounts.txt                                        ||
## - output will be counts in peaks called on sex chromosomes for each sample                      ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - shifted TagAlign files, split chromosomes for MACS3 TA in ${ALIGNED_DIR}/sexChr/              ||
## - MACS3 and BEDtools installed in a conda environment                                           ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR, XCHRBED,PEAKCOUNTS, PAR                      ||
##                                                                                                 ||
## ================================================================================================##

## ========================== ##
## MACS3 TA mode peak calling ##
## ========================== ##
  
for chr in chrX chrY
do

  echo "Starting peak calling using MACS3 TA for peaks in ${chr} at:"
  date -u
  f_TA=($(ls ${ALIGNED_DIR}/sexChr/*${chr}.tn5.tagAlign.gz))
  echo ${#f_TA[@]}
  cd ${PEAK_DIR}/ShiftedTagAlign/sexChr
  
  echo "Cutoff for broad peak calling is $MACS_CHR"
  macs3 callpeak -t ${f_TA[@]} -n ${chr} -f BED --outdir ${PEAK_DIR}/ShiftedTagAlign/sexChr -g 2.9e9 -q $MACS_CHR --keep-dup all --shift 100 --extsize 200 --nomodel --broad --broad-cutoff $MACS_CHR
   
  if [ ${chr} == chrX ]
	then 
		## X chr extract peaks nearest XIST and FIRRE
		bedtools closest -D a -id -io -a ${XCHRBED} -b ${PEAK_DIR}/ShiftedTagAlign/sexChr/chrX_peaks.broadPeak > ${PEAK_DIR}/ShiftedTagAlign/sexChr/chrX.broadPeak.filt		
	else
		## Y chr exclude peaks overlapping psuedoautosomal regions
	  bedtools intersect -v -a ${PEAK_DIR}/ShiftedTagAlign/sexChr/chrY_peaks.broadPeak -b ${PAR} \
    | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    | grep -P 'chr[\dXY]+[ \t]' > ${PEAK_DIR}/ShiftedTagAlign/sexChr/chrY.broadPeak.filt
	fi
  rm *_peaks.*
  
  ## count reads in peaks
	touch ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/${chr}.peakcounts.txt
	for file in ${f_TA[@]}
	do
		bedtools intersect -C -filenames -a ${PEAK_DIR}/ShiftedTagAlign/sexChr/${chr}.broadPeak.filt -b ${file} | awk -v var=${file} '{ print $0, var}' >> ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/${chr}.peakcounts.txt
    done
done

if [[ ! -f ${PEAK_DIR}/ShiftedTagAlign/sexChr/chrX.broadPeak.filt	]] || [[ ! -f ${PEAK_DIR}/ShiftedTagAlign/sexChr/chrY.broadPeak.filt ]] || [[ ! -f ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/chrX.peakcounts.txt ]] || [[ ! -f ${PEAKCOUNTS}/ShiftedTagAlign/sexChr/chrX.peakcounts.txt ]]
then
  { echo "Peaks in sex chromosomes were not successfully called. Please check STEP 5.2 SPLIT to make sure sex chromosomes reads exist." ; exit 1 ;}
else
  echo "Calling peaks in sex chromosomes using MACS3 TA completed."
  echo Job finished on:
  date -u
fi