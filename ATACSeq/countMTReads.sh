## get the number of reads aligned to MT chromosome from samtools idxstats output
## also count total number of mapped reads

echo "Filename\tMTReads\tAllMappedReads\n" > CountMTReads.txt

for file in ${ALIGNEDDIR}/*_statsperchr.txt
do 
    echo -n $(basename ${file}) >> CountMTReads.txt
    awk 'BEGIN{tot=0;mt=0} {tot=tot+$3} ($1 == "chrM") {mt=$3} END{printf "\t%d\t%d\n",mt,tot}' ${file} >> CountMTReads.txt
done


