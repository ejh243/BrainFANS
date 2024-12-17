#!/bin/bash
## =============================================================================================##
##                             ATAC-seq pipeline STEP 8.1: Collate QMetrics                     ##
## =============================================================================================##
## EXECUTION: sbatch ./subScripts/progressReportS2.sh                                           ||
## - execute from pipeline's main directory                                                     ||
##                                                                                              ||
## DESCRIPTION: This script  checks whether all samples have undergone all steps.               ||
##                                                                                              ||
## OUTPUT:                                                                                      ||
## - ${META_DIR}/summariseSampleProgressS2.csv                                                  ||
##                                                                                              ||
## REQUIRES:                                                                                    ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR,META_DIR                                   ||
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

## check for individual samples
## save output in csv file
echo "sampleID,sexChrSampleX,sexChrSampleY,GenotypeCheck" > ${META_DIR}/summariseSampleProgressS2.csv

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
    echo -n ${sampleName}, >> ${META_DIR}/summariseSampleProgressS2.csv

    if [ ! -s ${ALIGNED_DIR}/sexChr/${sampleName}.chrX.tn5.tagAlign.gz ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProgressS2.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProgressS2.csv
    fi
    
    if [ ! -s ${ALIGNED_DIR}/sexChr/${sampleName}.chrY.tn5.tagAlign.gz ]
    then
        echo -n "N," >> ${META_DIR}/summariseSampleProgressS2.csv
    else
        echo -n "Y," >> ${META_DIR}/summariseSampleProgressS2.csv
    fi
    
    if [ ! -s ${ALIGNED_DIR}/genotypeConcordance/${sampleName}*.log ]
    then
        echo "N" >> ${META_DIR}/summariseSampleProgressS2.csv
    else
        echo "Y" >> ${META_DIR}/summariseSampleProgressS2.csv
    fi
    

done

