#!/bin/bash
## ====================================================================================================##
##                   PIPPA STEP 4.2: Count MT chromosome reads and collate ENCODE results              ##
## ====================================================================================================##
## EXECUTION: sbatch ./subScripts/countMTReads.sh                                                      ||
##                                                                                                     ||
## DESCRIPTION: This script gets the number of reads aligned to MT chromosome and collates results     ||
## of ENCODE metrics calculation.                                                                      ||
##                                                                                                     ||
## OUTPUTS:                                                                                            ||
## ${ALIGNED_DIR}/countMTReads.txt                                                                     ||
## ${ALIGNED_DIR}/ENCODEMetrics/collateFlagStatMetrics.txt                                             ||
##                                                                                                     ||
## REQUIRES:                                                                                           ||
## - Variable in config file: ALIGNED_DIR                                                              ||
##                                                                                                     ||
## ====================================================================================================##


## ====================== ##
## COLLATE MT reads FILES ##
## ====================== ##

## Results will be output in txt file
echo -e "Filename\tMTReads\tAllMappedReads" > "${ALIGNED_DIR}/QCOutput/countMTReads.txt"

for file in ${ALIGNED_DIR}/*_statsperchr.txt
do 
    echo -n $(basename ${file%_statsperchr.txt}) >> "${ALIGNED_DIR}/QCOutput/countMTReads.txt"
    awk 'BEGIN{tot=0;mt=0} {tot=tot+$3} ($1 == "chrM") {mt=$3} END{printf "\t%d\t%d\n",mt,tot}' ${file} >> "${ALIGNED_DIR}/QCOutput/countMTReads.txt"
done

## ====================== ##
## COLLATE FLAGSTAT FILES ##
## ====================== ##

cd ${ALIGNED_DIR}/ENCODEMetrics/

for file in *flagstat.qc;
do
	echo -n ${file%.flagstat.qc},
   awk 'BEGIN { ORS = "," } {print $1} END { printf( "\n" ); }' ${file}
   
done  > "${ALIGNED_DIR}/QCOutput/collateFlagStatMetrics.txt"

## ====================== ##
## COLLATE ENCODE FILES   ##
## ====================== ##

cd ${ALIGNED_DIR}/ENCODEMetrics/

## Results will be output in txt file
echo -e "Filename\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2" > "${ALIGNED_DIR}/QCOutput/collateENCODEMetrics.txt"
for file in *pbc.qc;
do
  awk 'FNR==1{print FILENAME "\t" $0}' ${file}
done  >> "${ALIGNED_DIR}/QCOutput/collateENCODEMetrics.txt"


if [[ ! -f "${ALIGNED_DIR}/QCOutput/countMTReads.txt" ]] || [[ ! -f "${ALIGNED_DIR}/QCOutput/collateFlagStatMetrics.txt" ]]
then
  { echo "Results could not be collated" ; exit 1;}
else
  echo "Results collated successfully."
  echo Job ended on:
  date -u
fi

