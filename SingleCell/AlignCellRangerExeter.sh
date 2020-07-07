#!/bin/sh
#PBS -V #export all enviroment variables to the batch job
#PBS -q sq #submit to the serial queue
#PBS -l walltime=24:00:00 # maximum wall time for the job
#PBS -A Research_Project-MRC190311 #research project to submit under
#PBS -l nodes=1:ppn=1 #Number of processors 
#PBS -m e -M a.n.dahir@exeter.ac.uk # email me at job completion 

cd /gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/AlignedData/ExeterSequencing

module load CellRanger/3.1.0

cellranger count --id=10064_83550_sox10 \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/ftp1.sequencing.exeter.ac.uk/10064/10064_83550_sox10 \
--sample=10064_83550_sox10 \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=10064_83550_TN \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/ftp1.sequencing.exeter.ac.uk/10064/10064_83550_TN \
--sample=10064_83550_TN \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=10064_1991 \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/ftp1.sequencing.exeter.ac.uk/10064/10064_1991 \
--sample=10064_1991 \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=10064_1993 \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/ftp1.sequencing.exeter.ac.uk/10064/10064_1993 \
--sample=10064_1993 \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=10064_Hek293 \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/ftp1.sequencing.exeter.ac.uk/10064/10064_Hek293 \
--sample=10064_Hek293 \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA

cellranger count --id=10064_SHSY5Y \
--fastqs=/gpfs/mrc0/projects/Research_Project-MRC190311/SingleCell/ftp1.sequencing.exeter.ac.uk/10064/10064_SHSY5Y \
--sample=10064_SHSY5Y \
--transcriptome=/gpfs/mrc0/projects/Research_Project-MRC190311/References/premRNA/ftp1.sequencing.exeter.ac.uk/refdata-cellranger-GRCh38-3.0.0/pre-mRNA/GRCh38-3.0.0.pre-mRNA


