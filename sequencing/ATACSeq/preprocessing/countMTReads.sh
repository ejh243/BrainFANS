## ====================================================================================================##
##                    ATAC-seq pipeline STEP 4.2: Count MT chromosome reads                            ##
## ====================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/countMTReads.sh                                    ||
## - execute from scripts directory                                                                    ||
##                                                                                                     ||
## DESCRIPTION: This script gets the number of reads aligned to MT chromosome                          ||
##                                                                                                     ||
## OUTPUTS:                                                                                            ||
## countMTReads.txt                                                                                    ||
##                                                                                                     ||
## REQUIRES:                                                                                           ||
## - ALIGNEDDIR                                                                                        ||
##                                                                                                     ||
## ====================================================================================================##

## Results will be output in txt file
echo "Filename\tMTReads\tAllMappedReads\n" > ${ALIGNEDDIR}/countMTReads.txt

for file in ${ALIGNEDDIR}/*_statsperchr.txt
do 
    echo -n $(basename ${file%_statsperchr.txt}) >> ${ALIGNEDDIR}/countMTReads.txt
    awk 'BEGIN{tot=0;mt=0} {tot=tot+$3} ($1 == "chrM") {mt=$3} END{printf "\t%d\t%d\n",mt,tot}' ${file} >> ${ALIGNEDDIR}/countMTReads.txt
done