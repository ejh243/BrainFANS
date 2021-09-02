cd ${PEAKDIR}
mkdir -p BlacklistFiltered/MACS2
mkdir -p BlacklistFiltered/EPIC2

## process MACS2 called peaks
narrowPeak=($(ls MACS2/*_peaks.narrowPeak))

for peaks in ${narrowPeak[@]}
do
	bedtools intersect -a ${peaks} -b ${REFDIR}/ENCODE/Blacklist/hg38-blacklist.v2.bed -v > BlacklistFiltered/${peaks}
done


## process EPIC2 called peaks

epicPeak=($(ls EPIC2/*.txt))

for peaks in ${epicPeak[@]}
do
	bedtools intersect -a ${peaks} -b ${REFDIR}/ENCODE/Blacklist/hg38-blacklist.v2.bed -v > BlacklistFiltered/${peaks}
done