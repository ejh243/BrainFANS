#!/bin/bash
## =============================================================================================##
##                             ATAC-seq pipeline STEP 4.2: Collate QMetrics                     ##
## =============================================================================================##
## EXECUTION: sbatch ./subScripts/progressReport.sh                                             ||
## - execute from pipeline's main directory                                                     ||
##                                                                                              ||
## DESCRIPTION: This script  checks whether all samples have undergone all steps.               ||
##                                                                                              ||
## OUTPUT:                                                                                      ||
## - ${META_DIR}/summariseSampleProcessingProgress.csv                                          ||
##                                                                                              ||
## REQUIRES:                                                                                    ||
## - Variables in config file: RAWDATADIR, FASTQCDIR, ALIGNED_DIR, PEAK_DIR,META_DIR            ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                   ||
##                                                                                              ||
## =============================================================================================##

echo "Running progress report"

## create array of all fastq files
cd ${RAWDATADIR}

if test -f ${META_DIR}/samples.txt;
then 
    ## create an array from the file
    mapfile -t SAMPLEIDS < ${META_DIR}/samples.txt 
else
    exit 1
fi


### First count number of output files for each part of the pipeline

echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""

## check for FASTQC output
## NB two files per sample
echo "Number of fastqc reports found " $(ls ${FASTQCDIR}/*_fastqc.zip | wc -l)

## check for FASTP output
echo "Number of fastp reports found " $(ls ${TRIM_DIR}/qc/*.json | wc -l)

## check for bowtie output
echo "Number of Bowtie reports found " $(ls ${ALIGNED_DIR}/*.bowtie.log | wc -l)

echo "Number of filtered aligned files found " $(ls ${ALIGNED_DIR}/*.filt.nodup.bam | wc -l)

## check for output of ENCODE QC metrics calculation
## first delete empty files 

find ${ALIGNED_DIR}/ENCODEMetrics/ -size  0 -print -delete

echo "Number of ENCODE QC metric output files found " $(ls ${ALIGNED_DIR}/ENCODEMetrics/*.pbc.qc | wc -l)

## check for peak calling output

echo "Number of MACS3 peak files (paired end) found " $(ls ${PEAK_DIR}/MACS/BAMPE/*.broadPeak.filt | wc -l)


## check for individual samples
## save output in csv file
echo "sampleID,dataFolder,R1Filename,R2Filename,FASTQCR1,FASTQCR2,FASTP,BOWTIE,filteredAligned,ENCODEMetrics,MACS3PeaksPE" > ${META_DIR}/summariseSampleProcessingProgress.csv

for sampleName in ${SAMPLEIDS[@]}
do 
    toProcess=($(find ${RAWDATADIR} -maxdepth 1 -name ${sampleName}'*'))
  
    ## sort the toProcess array so that R1 and R2 are consecutive 
    IFS=$'\n' # need to set this as \n rather than default - a space, \t and then \n - so that elements are expanded using \n as delimiter
    toProcess=($(sort <<<"${toProcess[*]}")) ## sort so that the first element is R1
    unset IFS
    
    echo "Processing" ${sampleName}
    f1=$(basename ${toProcess[0]}) 
    f2=$(basename ${toProcess[1]})
    echo -n ${sampleName},${RAWDATADIR}, ${f1},${f2}, >> ${META_DIR}/summariseSampleProcessingProgress.csv

    if [ ! -s ${FASTQCDIR}/${f1%%.*}*fastqc.zip ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${FASTQCDIR}/${f2%%.*}*fastqc.zip ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi

    if [ ! -s ${TRIM_DIR}/qc/${sampleName}*.json ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ $(wc -l < ${ALIGNED_DIR}/${sampleName}.bowtie.log) != 15 ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNED_DIR}/${sampleName}.filt.nodup.bam ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNED_DIR}/ENCODEMetrics/${sampleName}*.pbc.qc ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${PEAK_DIR}/MACS/BAMPE/${sampleName}.broadPeak.filt ]
    then
        echo "N" >> ${META_DIR}/summariseSampleProcessingProgress.csv
    else
        echo "Y" >> ${META_DIR}/summariseSampleProcessingProgress.csv
    fi
    

done

