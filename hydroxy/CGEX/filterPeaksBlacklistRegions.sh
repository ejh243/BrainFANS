cd ${PEAKDIR}

mkdir -p BlacklistFiltered

narrowPeak=($(ls *_peaks.narrowPeak))

for peaks in ${narrowPeak[@]}
do
	bedtools intersect -a ${peaks} -b ${REFDIR}/ENCODE/Blacklist/hg38-blacklist.v2.bed -v > BlacklistFiltered/${peaks}
done