#!/bin/bash
## ================================================================================================##
##                                   PIPPA STEP 5.2: Sex Chr Peak calling                          ||                  
## ================================================================================================##
## EXECUTION: sbatch ./subScripts/sexChrPeaks.sh                                                   ||
##                                                                                                 ||
## DESCRIPTION: This script call peaks only on reads assigned to sex chromosomes merging all       ||
## samples, using MACS3 paired-end mode                                                            ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## PEAK_DIR/chrX.narrowPeak.filt, PEAK_DIR/chrY.narrowPeak.filt                                    ||
## - output will be peaks called in sex chromosomes and not aligned to PAR regions                 ||
##  PEAKCOUNTS/chrX.peakcounts.txt, PEAKCOUNTS/chrY.peakcounts.txt                                 ||
## - output will be counts in peaks called on sex chromosomes for each sample                      ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - bam files, split chromosomes for MACS3 in ALIGNED_DIR/sexChr/                                 ||
## - MACS3, bedtools installed in conda environment                                                ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR, XCHRBED,PEAKCOUNTS, PAR                      ||
##                                                                                                 ||
## ================================================================================================##


## ========================== ##
## MACS3 PE mode peak calling ##
## ========================== ##

for chr in chrX chrY
do
  echo "Calling peaks in ${chr} using MACS3 PE"
  mapfile -t SAMPLES < ${META_DIR}/samples.txt 
  
  f_PE=($(ls ${ALIGNED_DIR}/sexChr/*.${chr}.bam))
  
  echo "Number of BAM files found for chromosome ${chr}: ${#f_PE[@]}"
  cd "${PEAK_DIR_SEXCHR}"
  macs3 callpeak -t ${f_PE[@]} -n ${chr} -f BAMPE --outdir "${PEAK_DIR_SEXCHR}" \
   -g 2.9e9 -q "${macs_sexchr_thres}" \
   --keep-dup all --nomodel 

  
	if [ ${chr} == chrX ]
	then 
		## X chr extract peaks nearest XIST and FIRRE
		bedtools closest -D a -id -io -a ${XCHRBED} -b "${PEAK_DIR_SEXCHR}/chrX_peaks.narrowPeak" > "${PEAK_DIR_SEXCHR}/chrX.narrowPeak.filt"		
	else
		## Y chr exclude peaks overlapping psuedoautosomal regions
	  bedtools intersect -v -a "${PEAK_DIR_SEXCHR}/chrY_peaks.narrowPeak" -b ${PAR} \
    	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    	| grep -P 'chr[\dXY]+[ \t]' > "${PEAK_DIR_SEXCHR}/chrY.narrowPeak.filt"
	fi
  rm *_peaks.*
 
  ## count reads in peaks
  
  if [[ -f "${COUNTS_DIR_SEXCHR}/${chr}.peakcounts.txt"	]]
  then
    rm "${COUNTS_DIR_SEXCHR}/${chr}.peakcounts.txt"
  fi
  
	touch "${COUNTS_DIR_SEXCHR}/${chr}.peakcounts.txt"
	for file in ${f_PE[@]}
	do
		bedtools intersect -C -filenames -a "${PEAK_DIR_SEXCHR}/${chr}.narrowPeak.filt" -b ${file} \
		| awk -v var=${file} '{ print $0, var}' >> "${COUNTS_DIR_SEXCHR}/${chr}.peakcounts.txt"
    done
done