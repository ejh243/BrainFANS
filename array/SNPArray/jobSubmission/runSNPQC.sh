#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=SNPQC_%j.log
#SBATCH --error=SNPQC_%j.err
#SBATCH --job-name=SNPQC


## print start date and time
echo Job started on:
date -u


source $1 || exit

echo "Processing data located in :" ${DATADIR}

cd ${SCRIPTDIR}/array/SNPArray/preprocessing

sh 1_QC.sh

module load $RVERS
sh 2_CheckEthnicity.sh


## create sample sheets of samples classed as different populations
## run relatedness check within ethnicites
populations=($(cut -f3 --delim="," ${PROCESSDIR}/PredictedPopulations.csv | tail -n +2 | sort | uniq))
for each in ${populations[@]}
do
   grep ${each} ${PROCESSDIR}/PredictedPopulations.csv | cut -f1-2 --delim="," --output-delimiter=" " > ${PROCESSDIR}/${each}Samples.txt
   
   ## 
   if [[ $(wc -l <${PROCESSDIR}/${each}Samples.txt) -ge 1 ]]
   then
      sh 3_CheckRelatedness.sh ${PROCESSDIR}/${each}Samples.txt ${each}
   fi
   
done   

## plot relatedness statistics acroos whole cohort 

Rscript utilitys/plotKinshipCoeff.r ${PROCESSDIR}/QCoutput/

module purge
module load VCFtools
sh 4_formatForImputation.sh ALL ${KGG}/1000GP_Phase3_combined.legend
sh 4_formatForImputation.sh EUR ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

## print finish date and time
echo Job finished on:
date -u


mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/SNPQC_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/SNPQC_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/SNPQC_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/SNPQC_${SLURM_JOB_ID}.err"