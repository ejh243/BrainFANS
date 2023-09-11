## parameter choices guided by (https://baderlab.github.io/CBW_Pathways_2021/lectures/Pathways_2021_segway_semi_automated_genome_annotation_post_submission_draft.pdf)

label=$1
shift
FILES=$@

echo
echo 'Starting annotating model at: '
date -u 

echo 'Input samples are:' $FILES

## annotate the whole genome using the trained model
segway annotate --bigBed=${ANNODIR}/segway.layered.bb ${FILES} \
	${TRAINDIR} ${ANNODIR}

if [[ $? == 0 ]]
then
	echo 'Model annotated'
	date -u
fi

## calculate and plot the length distribution of segments in each label, and the genomic fraction covered by each label
cd ${ANNODIR}
segtools-length-distribution segway.bed.gz


## calculate the enrichment of each segment label over a gene annotation 
version=${GENCODE%.annotation.gtf.gz}
echo 'Calculating enrichment using gene annotation:' $(basename ${version})

# This command runs segtools-aggregation in “gene” mode and creates two figures showing
# the enrichment of a segmentation over an idealized transcriptional (feature_aggregation.splicing.png)
# and translational (feature_aggregation.translation.png) gene model.

segtools-aggregation --mode=gene --normalize --outdir ${AGGREDIR} \
	${ANNODIR}/segway.bed.gz ${GENCODE}

## Create a filtered GENCODE annotation list with only genes if it doesn't already exist
if ! test -f ${version}.genes.gtf
then
	zcat ${GENCODE} | sed '/\(gene\t\|^#\)/!d' > ${version}.genes.gtf
fi

let "label -= 1"
echo $label
## Create a list of GENCODE genes that overlap with each label from the produced annotation
for ((i=0;i<=label;i++));
do 
	zcat ${ANNODIR}/segway.bed.gz | awk --assign label=${i} 'BEGIN{OFS="\t"} { if ($4==label) print $1,$2,$3 }' | bedtools \
		intersect -a ${version}.genes.gtf -b stdin | sed 's/.*gene_id "\([^"]\+\).*/\1/' | uniq > "${AGGREDIR}/segway.${i}.$(basename ${version}).genes.bed"
done

if [[ $? == 0 ]]
then
	echo 'Gene enrichment aggregated '
	date -u
fi