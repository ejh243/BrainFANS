## parameter choices guided by (https://baderlab.github.io/CBW_Pathways_2021/lectures/Pathways_2021_segway_semi_automated_genome_annotation_post_submission_draft.pdf)
state=$1
shift
FILES=$@

echo $state

echo $FILES

echo
echo "Starting training Segway at: "
date -u 

echo "Input samples are:" $FILES
echo "traindir =" ${TRAINDIR}

# train segway on 1% of the genome, creating a 10 label model, using 10 simultaneous training instances 
segway train --resolution=10 --num-instances=10 --minibatch-fraction=0.01 --num-labels=${state} \
	--max-train-rounds=10 ${FILES} ${TRAINDIR}

cd $TRAINDIR

# Plot the emission parameters learned during the training task
segtools-gmtk-parameters ${TRAINDIR}/params/params.params


if [[ $? == 0 ]]
then
	echo "Model trained"
	date -u
fi