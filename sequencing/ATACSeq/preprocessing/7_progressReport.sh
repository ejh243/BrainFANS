## searches for output files from each part of the QC pipeline and reports number of files found
## then searches specifically for each sample to identify what is missing

## EXECUTION
# sh ./ATACSeq/preprocessing/7_progressReport.sh
# assumes config file has been loaded
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# RAWDATADIR, FASTQCDIR, ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# -

## INPUT
# -

## OUTPUT
# ${METADIR}/summariseSampleProcessingProgress.csv


echo "Running progress report"

## create array of all fastq files
cd ${RAWDATADIR}

if test -f ${METADIR}/samples.txt;
then 
    ## create an array from the file
    mapfile -t SAMPLEIDS < ${METADIR}/samples.txt 
else
    exit 1
fi


### First count number of output files for each part of the pipeline

echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""

## check for FASTQC output
## NB two files per sample
echo "Number of fastqc reports found " $(ls ${FASTQCDIR}/*_fastqc.zip | wc -l)

## check for FASTP output
echo "Number of fastp reports found " $(ls ${TRIMDIR}/fastp_reports/*.json | wc -l)

## check for bowtie output
echo "Number of Bowtie reports found " $(ls ${ALIGNEDDIR}/*.bowtie.log | wc -l)

echo "Number of filtered aligned files found " $(ls ${ALIGNEDDIR}/*.filt.nodup.bam | wc -l)

## check for output of ENCODE QC metrics calculation
## first delete empty files 

find ${ALIGNEDDIR}/ENCODEMetrics/ -size  0 -print -delete

echo "Number of ENCODE QC metric output files found " $(ls ${ALIGNEDDIR}/ENCODEMetrics/*.pbc.qc | wc -l)


## check for peak calling output

echo "Number of MACS2 peak files (shifted tag align) found " $(ls ${PEAKDIR}/MACS/ShiftedTagAlign/*.broadPeak.filt | wc -l)

echo "Number of MACS2 peak files (paired end) found " $(ls ${PEAKDIR}/MACS/BAMPE/*.broadPeak.filt | wc -l)


## check for individual samples
## save output in txt file
echo "sampleID,dataFolder,R1Filename,R2Filename,FASTQCR1,FASTQCR2,FASTP,BOWTIE,filteredAligned,ENCODEMetrics,MACS2Peaks1,MACS2Peaks2" > ${METADIR}/summariseSampleProcessingProgress.csv

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

    echo -n ${sampleName},${RAWDATADIR}, ${f1},${f2}, >> ${METADIR}/summariseSampleProcessingProgress.csv

    
    if [ ! -s ${FASTQCDIR}/${f1%%.*}*fastqc.zip ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${FASTQCDIR}/${f2%%.*}*fastqc.zip ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi

    if [ ! -s ${TRIMDIR}/fastp_reports/${sampleName}*.json ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNEDDIR}/${sampleName}.bowtie.log ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNEDDIR}/ENCODEMetrics/${sampleName}*.pbc.qc ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt ]
    then
        echo -n "N," >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt ]
    then
        echo "N" >> ${METADIR}/summariseSampleProcessingProgress.csv
    else
        echo "Y" >> ${METADIR}/summariseSampleProcessingProgress.csv
    fi
  
done
