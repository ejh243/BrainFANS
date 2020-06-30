#!/bin/sh
#PBS -V #export all enviroment variables to the batch job
#PBS -q sq #submit to the serial queue
#PBS -l walltime=24:00:00 # maximum wall time for the job
#PBS -A Research_Project-MRC190311 #research project to submit under
#PBS -l nodes=1:ppn=1 #Number of processors 
#PBS -m e -M a.n.dahir@exeter.ac.uk # email me at job completion 

cd /gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/AlignedData/

module load CellRanger/3.1.0

cellranger count --id=AN17142-TN \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/RawData/SourceBiosciences/Jonathan_Mill_SOUK005843/AN17142-TN \
--sample=AN17142-TN \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=AN17142-Sox10 \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/RawData/SourceBiosciences/Jonathan_Mill_SOUK005843/AN17142-Sox10 \
--sample=AN17142-Sox10 \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=AN17142-NeuN \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/RawData/SourceBiosciences/Jonathan_Mill_SOUK005843/AN17142-NeuN \
--sample=AN17142-NeuN \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA