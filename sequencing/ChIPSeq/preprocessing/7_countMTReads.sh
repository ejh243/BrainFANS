## get the number of reads aligned to MT chromosome from samtools idxstats output


## EXECUTION
# sh ATACSeq/countMTReads.sh 

## REQUIRES the following variables in config file
# ALIGNEDDIR

## REQUIRES the following software
# -

## INPUT
# -

## OUTPUT
# ${ALIGNEDDIR}/countMTReads.txt


echo "Filename\tMTReads\tAllMappedReads\n" > ${ALIGNEDDIR}/countMTReads.txt

for file in ${ALIGNEDDIR}/*_statsperchr.txt
do 
    echo -n $(basename ${file%_statsperchr.txt}) >> ${ALIGNEDDIR}/countMTReads.txt
    awk 'BEGIN{tot=0;mt=0} {tot=tot+$3} ($1 == "chrM") {mt=$3} END{printf "\t%d\t%d\n",mt,tot}' ${file} >> ${ALIGNEDDIR}/countMTReads.txt
done