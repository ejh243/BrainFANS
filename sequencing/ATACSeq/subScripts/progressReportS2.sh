#!/bin/bash
## =============================================================================================##
##                           PIPPA STEP 7.1: Progress Report Stage 2 QC                         ##
## =============================================================================================##
## EXECUTION: sbatch ./subScripts/progressReportS2.sh                                           ||
##                                                                                              ||
## DESCRIPTION: This script checks whether all samples have undergone all steps in Stage 2 QC.  ||
##                                                                                              ||
## OUTPUT:                                                                                      ||
## - ${META_DIR}/summariseSampleProgressS2.csv                                                  ||
##                                                                                              ||
## REQUIRES:                                                                                    ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR,META_DIR                                   ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                   ||
##                                                                                              ||
## =============================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

echo "Running progress report"

## create array of all fastq files
cd ${RAWDATADIR}

if test -f "${SAMPLES_LIST}";
then 
    ## create an array from the file
    mapfile -t SAMPLEIDS < "${SAMPLES_LIST}"
else
    exit 1
fi

### First count number of output files for each part of the pipeline

echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""

## ============ ##
##    COLLATE   ##
## ============ ##

## check for individual samples
## save output in csv file
echo "sampleID,sexChrSampleX,sexChrSampleY,GenotypeCheck" > "${META_DIR}/summariseSampleProgressS2.csv"

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
    echo -n ${sampleName}, >> "${META_DIR}/summariseSampleProgressS2.csv"

    if [ ! -s "${ALIGNED_DIR}/sexChr/${sampleName}.chrX.bam" ]
    then
        echo -n "N," >> "${META_DIR}/summariseSampleProgressS2.csv"
    else
        echo -n "Y," >> "${META_DIR}/summariseSampleProgressS2.csv"
    fi
    
    if [ ! -s "${ALIGNED_DIR}/sexChr/${sampleName}.chrY.bam" ]
    then
        echo -n "N," >> "${META_DIR}/summariseSampleProgressS2.csv"
    else
        echo -n "Y," >> "${META_DIR}/summariseSampleProgressS2.csv"
    fi
    
    if [ ! -s ${ALIGNED_DIR}/genotypeConcordance/${sampleName}*.log ]
    then
        echo "N" >> "${META_DIR}/summariseSampleProgressS2.csv"
    else
        echo "Y" >> "${META_DIR}/summariseSampleProgressS2.csv"
    fi

done

if [[ ! -f "${META_DIR}/summariseSampleProgressS2.csv" ]]
then
  { echo "Stage 2 progress report has been not been created." ; exit 1; }
else
  echo "Stage 2 progress report has been been created."
  echo Job ended on:
  date -u
fi