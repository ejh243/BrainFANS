toProcess=$1
sampleName=${toProcess%.bam}

module load BEDTools

# Create virtual SE file containing both read pairs
bedtools bamtobed -i $toProcess | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${sampleName}.PE2SE.tn5.tagAlign.gz

# ================
# Shift tagAlign file
# ================

zcat ${sampleName}.PE2SE.tn5.tagAlign.gz | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${sampleName}.tn5.tagAlign.gz

rm ${sampleName}.PE2SE.tn5.tagAlign.gz