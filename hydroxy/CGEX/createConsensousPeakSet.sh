## comparison of methods to call peaks found that some cell-type spefic peaks were missed when calling across all samples, but calling across all samples gained some lower level peaks.
## hence we will combine peaks called within each cell type and brain region and those call across all to create the consensous peak set
## peaks have already been filtered to exclude those overlapping the blacklist regions
## need to filter peaks based on q-value first.
## do this for all higher level categorisations, e.g. brain regions and cell types


cd ${PEAKDIR}/BlacklistFiltered

mkdir -p ConsensousPeakSets

## create list of categories
groups=($(ls | grep '^[[:upper:]]*[[:digit:]]*_peaks.narrowPeak'))

### identify cell type specific peaks (i.e. those unique to a cell type and not in all peak set)

for all in ${groups[@]}
do
  label=${all%_peaks.narrowPeak}
  ctPeaks=($(ls *${label}*_peaks.narrowPeak))
  for type in ${ctPeaks[@]}
  do
    if [ ${type} != ${all} ]
    then
      bedtools intersect -a ${type} -b ${all} -v > unique_${type}
    fi
  done

  ## merge with peaks called across all samples
 
  uniqueCTPeaks=($(ls unique_*))
  cat ${all} ${uniqueCTPeaks[@]} | bedtools -i sort > ConsensousPeakSets/merged_${all}

done


