## comparison of methods to call peaks found that some cell-type spefic peaks were missed when calling across all samples, but calling across all samples gained some lower level peaks.
## hence we will combine peaks called within each cell type and brain region and those call across all to create the consensous peak set
## peaks have already been filtered to exclude those overlapping the blacklist regions
## need to filter peaks based on q-value first.
## do this for all higher level categorisations, e.g. brain regions and cell types




mkdir -p ${PEAKDIR}/ConsensousPeakSets/MACS2
mkdir -p ${PEAKDIR}/ConsensousPeakSets/EPIC2

## MACS2 called peaks

cd ${PEAKDIR}/BlacklistFiltered/MACS2/

## create list of categories
groups=($(ls | grep '^[[:upper:]]*[[:digit:]]*_peaks.narrowPeak'))

### identify cell type specific peaks (i.e. those unique to a cell type and not in all peak set)

for all in ${groups[@]}
do
  label=${all%_peaks.narrowPeak}
  ctPeaks=($(ls *${label}*_peaks.narrowPeak))
  
  ## filter peaks by q-value
  awk '{if ($9 > 1.30) print $0}' $all > tmp_$all
  for type in ${ctPeaks[@]}
  do
    if [ ${type} != ${all} ]
    then
	  awk '{if ($9 > 1.30) print $0}' $type > tmp_$type
      bedtools intersect -a tmp_${type} -b tmp_${all} -v > unique_${type}
    fi
  done

  ## merge with peaks called across all samples
 
  uniqueCTPeaks=($(ls unique_*))
  cat tmp_${all} ${uniqueCTPeaks[@]} | bedtools sort -i > ${PEAKDIR}/ConsensousPeakSets/MACS2/merged_${all}
  rm tmp_*
  rm unique_*
done

## EPIC2 called peaks
cd ${PEAKDIR}/BlacklistFiltered/EPIC2/

## create list of categories
groups=($(ls | grep '^[[:upper:]]*[[:digit:]]*.txt'))

### identify cell type specific peaks (i.e. those unique to a cell type and not in all peak set)

for all in ${groups[@]}
do
  label=${all%.txt}
  ctPeaks=($(ls *${label}.txt))
  
  ## filter peaks by q-value
  awk '{if ($9 < 0.05) print $0}' $all > tmp_$all
  for type in ${ctPeaks[@]}
  do
    if [ ${type} != ${all} ]
    then
	  awk '{if ($9 < 0.05) print $0}' $type > tmp_$type
      bedtools intersect -a tmp_${type} -b tmp_${all} -v > unique_${type}
    fi
  done

  ## merge with peaks called across all samples
 
  uniqueCTPeaks=($(ls unique_*))
  cat tmp_${all} ${uniqueCTPeaks[@]} | bedtools sort -i > ${PEAKDIR}/ConsensousPeakSets/EPIC2/merged_${all}
  rm tmp_*
  rm unique_*
done
