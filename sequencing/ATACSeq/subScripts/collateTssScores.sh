#!/bin/bash
## ===================================================================================================================##
##                                   PIPPA STEP 2.4: TSS enrichment results                                           ##
## ===================================================================================================================##
## EXECUTION: sbatch ./subScripts/collateTSSscores.sh                                                                 ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script collates the results of TSS enrichment scores (STEP 2.2-2.3) single csv file.             ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - "ALIGNED_DIR/QCOUTPUT/TSS_enrichment.tsv"                                                                        ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Variables in config file: ALIGNED_DIR                                                                            ||
## -_TSS.mat.gz TSS enrichment score files for all samples                                                            ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

TSSFILE="${ALIGNED_DIR}/QCOutput/TSS_enrichment.tsv"

## ============ ##
##    COLLATE   ##
## ============ ##

echo -e "sample\tTSS_enrichment" > "${TSSFILE}"

cd ${TSS_DIR}
for mat in ${TSS_DIR}/*_TSS.mat.gz; do
  sample=$(basename "$mat" _TSS.mat.gz)
  zcat "$mat" | awk -v sample="$sample" '
  BEGIN { center=0; flank=0; n=0 }
    NR>1 {
      c=0; f=0
      # TSS center ±50 bp (bins 195–205 ? cols 201–211)
      for (i=201;i<=211;i++) c += $i
      # flanks (bins 2–101 and 301–400 ? cols 8–107 and 307–406)
      for (i=8;i<=107;i++) f += $i
      for (i=307;i<=406;i++) f += $i
      center+= c/11
      flank+= f/200
      n++
    }
    END {
    print sample "\t" center/flank
    }
  ' >> "${TSSFILE}"
done

if [[ ! -f "${TSSFILE}" ]]
then
  { echo "TSS enrichment scores results could not be collated. Please check STEP 2.2-2.3." ; exit 1;}
else
  echo "TSS enrichment scores results have been collated."
  echo Job finished on:
  date -u
fi