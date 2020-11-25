
cd ${ALIGNEDDIR}

## check if required output directories exist, if not create
outdir=$PEAKDIR/macs2.narrow.aug18 #for macs2
outdir2=$PEAKDIR/macs2.narrow.aug18.dedup #for macs2 dedup version

outdirbroad=$PEAKDIR/macs2.broad.aug18 #for macs2
outdirbroad2=$PEAKDIR/macs2.broad.aug18.dedup #for macs2 dedup version

outdirseac=$PEAKDIR/seacr.aug12 #for seacr
outdirseac2=$PEAKDIR/seacr.aug12.dedup #for seacr dedup version

for d in $outdir $outdir2 $outdirbroad $outdirbroad2 $outdirseac $outdirseac2; do
if [ ! -d $d ]; then
	mkdir $d
fi
done

outdir=$PEAKDIR/macs2.narrow.all.frag.aug18 #for macs2
outdir2=$PEAKDIR/macs2.narrow.all.frag.aug18.dedup #for macs2 dedup version

outdirbroad=$PEAKDIR/macs2.broad.all.frag.aug18 #for macs2
outdirbroad2=$PEAKDIR/macs2.broad.all.frag.aug18.dedup #for macs2 dedup version

#SEACR peak calling
outdirseac=$PEAKDIR/seacr.aug12.all.frag #for seacr
outdirseac2=$PEAKDIR/seacr.aug12.all.frag.dedup #for seacr dedup version

for d in $outdir $outdir2 $outdirbroad $outdirbroad2 $outdirseac $outdirseac2; do
if [ ! -d $d ]; then
mkdir $d
fi
done



BAMFILES=($(ls *_aligned_reads.bam))
for f in ${BAMFILES[@]}
do
	sampleName=${f/_aligned_reads.bam}
		
	echo "Peak calling using MACS2... ""$sampleName".bam
	date
	bam_file=dup.marked.120bp/"$sampleName".bam
	outdir=$PEAKDIR/macs2.narrow.aug18 #for macs2
	outdir2=$PEAKDIR/macs2.narrow.aug18.dedup #for macs2 dedup version

	outdirbroad=$PEAKDIR/macs2.broad.aug18 #for macs2
	outdirbroad2=$PEAKDIR/macs2.broad.aug18.dedup #for macs2 dedup version

	outdirseac=$PEAKDIR/seacr.aug12 #for seacr
	outdirseac2=$PEAKDIR/seacr.aug12.dedup #for seacr dedup version

	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $outdir/"$sampleName".macs2
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdir2 -q 0.01 -B --SPMR 2> $outdir2/"$sampleName".dedup.macs2

	#broad peak calls
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> $outdirbroad/"$sampleName".broad.all.frag.macs2
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirbroad2 --broad --broad-cutoff 0.1 -B --SPMR 2> $outdirbroad2/"$sampleName".broad.all.frag.dedup.macs2
	
	python $CUTRUNPATH/get_summits_broadPeak.py $outdirbroad/"$sampleName"_peaks.broadPeak|sort-bed - > $outdirbroad/"$sampleName"_summits.bed
	python $CUTRUNPATH/get_summits_broadPeak.py $outdirbroad2/"$sampleName"_peaks.broadPeak|sort-bed - > $outdirbroad2/"$sampleName"_summits.bed

	#SEACR peak calls
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirseac -q 0.01 -B --keep-dup all
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirseac2 -q 0.01 -B
	python $CUTRUNPATH/change.bdg.py $outdirseac/"$sampleName"_treat_pileup.bdg > $outdirseac/"$sampleName"_treat_integer.bdg
	python $CUTRUNPATH/change.bdg.py $outdirseac2/"$sampleName"_treat_pileup.bdg > $outdirseac2/"$sampleName"_treat_integer.bdg
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac/"$sampleName"_treat_integer.bdg 0.01 non stringent $outdirseac/"$sampleName"_treat $Rscriptbin
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac2/"$sampleName"_treat_integer.bdg 0.01 non stringent $outdirseac2/"$sampleName"_treat $Rscriptbin
	sort-bed $outdirseac/"$sampleName"_treat.stringent.bed > $outdirseac/"$sampleName"_treat.stringent.sort.bed
	sort-bed $outdirseac2/"$sampleName"_treat.stringent.bed > $outdirseac2/"$sampleName"_treat.stringent.sort.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac/"$sampleName"_treat.stringent.bed|sort-bed - > $outdirseac/"$sampleName"_treat.stringent.sort.summits.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac2/"$sampleName"_treat.stringent.bed|sort-bed - > $outdirseac2/"$sampleName"_treat.stringent.sort.summits.bed
	for i in _summits.bed _peaks.xls _peaks.narrowPeak _control_lambda.bdg _treat_pileup.bdg; do 
	rm -rf $outdirseac/"$sampleName"$i
	rm -rf $outdirseac2/"$sampleName"$i
	done


	#SEACR relaxed peak calls
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac/"$sampleName"_treat_integer.bdg 0.01 non relaxed $outdirseac/"$sampleName"_treat $Rscriptbin
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac2/"$sampleName"_treat_integer.bdg 0.01 non relaxed $outdirseac2/"$sampleName"_treat $Rscriptbin
	sort-bed $outdirseac/"$sampleName"_treat.relaxed.bed > $outdirseac/"$sampleName"_treat.relaxed.sort.bed
	sort-bed $outdirseac2/"$sampleName"_treat.relaxed.bed > $outdirseac2/"$sampleName"_treat.relaxed.sort.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac/"$sampleName"_treat.relaxed.bed|sort-bed - > $outdirseac/"$sampleName"_treat.relaxed.sort.summits.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac2/"$sampleName"_treat.relaxed.bed|sort-bed - > $outdirseac2/"$sampleName"_treat.relaxed.sort.summits.bed

	cur=`pwd`
	echo "Converting bedgraph to bigwig... ""$sampleName".bam
	date
#	cd $outdir
	LC_ALL=C sort -k1,1 -k2,2n $outdir/"$sampleName"_treat_pileup.bdg > $outdir/"$sampleName".sort.bdg
	## need to exclude chrEBV entries as error otherwise
	awk '{if ($1 != "chrEBV") print $0}' $outdir/"$sampleName".sort.bdg > $outdir/"$sampleName".tmp.sort.bdg
	$CUTRUNPATH/bedGraphToBigWig $outdir/"$sampleName".tmp.sort.bdg $chromsizedir/hg38.chrom.sizes $outdir/"$sampleName".sorted.bw
	rm -rf $outdir/"$sampleName".sort.bdg
	rm -rf $outdir/"$sampleName".tmp.sort.bdg

	#cd $outdir2
	LC_ALL=C sort -k1,1 -k2,2n $outdir2/"$sampleName"_treat_pileup.bdg > $outdir2/"$sampleName".sort.bdg
	awk '{if ($1 != "chrEBV") print $0}' $outdir2/"$sampleName".sort.bdg > $outdir2/"$sampleName".tmp.sort.bdg
	$CUTRUNPATH/bedGraphToBigWig $outdir2/"$sampleName".tmp.sort.bdg $chromsizedir/hg38.chrom.sizes $outdir2/"$sampleName".sorted.bw
	rm -rf $outdir2/"$sampleName".sort.bdg
	rm -rf $outdir2/"$sampleName".tmp.sort.bdg

	#all fragments
	bam_file=dup.marked/"$sampleName".bam

	outdir=$PEAKDIR/macs2.narrow.all.frag.aug18 #for macs2
	outdir2=$PEAKDIR/macs2.narrow.all.frag.aug18.dedup #for macs2 dedup version

	outdirbroad=$PEAKDIR/macs2.broad.all.frag.aug18 #for macs2
	outdirbroad2=$PEAKDIR/macs2.broad.all.frag.aug18.dedup #for macs2 dedup version

	#SEACR peak calling
	outdirseac=$PEAKDIR/seacr.aug12.all.frag #for seacr
	outdirseac2=$PEAKDIR/seacr.aug12.all.frag.dedup #for seacr dedup version


	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $outdir/"$sampleName".macs2
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdir2 -q 0.01 -B --SPMR 2> $outdir2/"$sampleName".dedup.macs2

	#broad peak calls
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --keep-dup all 2> $outdirbroad/"$sampleName".broad.all.frag.macs2
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirbroad2 --broad --broad-cutoff 0.1 -B 2> $outdirbroad2/"$sampleName".broad.all.frag.dedup.macs2
	python $CUTRUNPATH/get_summits_broadPeak.py $outdirbroad/"$sampleName"_peaks.broadPeak|sort-bed - > $outdirbroad/"$sampleName"_summits.bed
	python $CUTRUNPATH/get_summits_broadPeak.py $outdirbroad2/"$sampleName"_peaks.broadPeak|sort-bed - > $outdirbroad2/"$sampleName"_summits.bed

	#SEACR peak calling
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirseac -q 0.01 -B --keep-dup all
	macs2 callpeak -t $bam_file -g hs -f BAMPE -n $sampleName --outdir $outdirseac2 -q 0.01 -B 
	python $CUTRUNPATH/change.bdg.py $outdirseac/"$sampleName"_treat_pileup.bdg > $outdirseac/"$sampleName"_treat_integer.bdg
	python $CUTRUNPATH/change.bdg.py $outdirseac2/"$sampleName"_treat_pileup.bdg > $outdirseac2/"$sampleName"_treat_integer.bdg
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac/"$sampleName"_treat_integer.bdg 0.01 non stringent $outdirseac/"$sampleName"_treat $Rscriptbin
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac2/"$sampleName"_treat_integer.bdg 0.01 non stringent $outdirseac2/"$sampleName"_treat $Rscriptbin
	sort-bed $outdirseac/"$sampleName"_treat.stringent.bed > $outdirseac/"$sampleName"_treat.stringent.sort.bed
	sort-bed $outdirseac2/"$sampleName"_treat.stringent.bed > $outdirseac2/"$sampleName"_treat.stringent.sort.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac/"$sampleName"_treat.stringent.bed|sort-bed - > $outdirseac/"$sampleName"_treat.stringent.sort.summits.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac2/"$sampleName"_treat.stringent.bed|sort-bed - > $outdirseac2/"$sampleName"_treat.stringent.sort.summits.bed
	for i in _summits.bed _peaks.xls _peaks.narrowPeak _control_lambda.bdg _treat_pileup.bdg; do 
	rm -rf $outdirseac/"$sampleName"$i
	rm -rf $outdirseac2/"$sampleName"$i
	done

	$CUTRUNPATH/SEACR_1.1.sh $outdirseac/"$sampleName"_treat_integer.bdg 0.01 non relaxed $outdirseac/"$sampleName"_treat $Rscriptbin
	$CUTRUNPATH/SEACR_1.1.sh $outdirseac2/"$sampleName"_treat_integer.bdg 0.01 non relaxed $outdirseac2/"$sampleName"_treat $Rscriptbin
	sort-bed $outdirseac/"$sampleName"_treat.relaxed.bed > $outdirseac/"$sampleName"_treat.relaxed.sort.bed
	sort-bed $outdirseac2/"$sampleName"_treat.relaxed.bed > $outdirseac2/"$sampleName"_treat.relaxed.sort.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac/"$sampleName"_treat.relaxed.bed|sort-bed - > $outdirseac/"$sampleName"_treat.relaxed.sort.summits.bed
	python $CUTRUNPATH/get_summits_seacr.py $outdirseac2/"$sampleName"_treat.relaxed.bed|sort-bed - > $outdirseac2/"$sampleName"_treat.relaxed.sort.summits.bed

	echo "Converting bedgraph to bigwig... ""$sampleName".bam
	date
	LC_ALL=C sort -k1,1 -k2,2n $outdir/"$sampleName"_treat_pileup.bdg > $outdir/"$sampleName".sort.bdg
	awk '{if ($1 != "chrEBV") print $0}' $outdir/"$sampleName".sort.bdg > $outdir/"$sampleName".tmp.sort.bdg
	$CUTRUNPATH/bedGraphToBigWig $outdir/"$sampleName".tmp.sort.bdg $chromsizedir/hg38.chrom.sizes $outdir/"$sampleName".sorted.bw
	rm -rf $outdir/"$sampleName".sort.bdg
	rm -rf $outdir/"$sampleName".tmp.sort.bdg


	LC_ALL=C sort -k1,1 -k2,2n $outdir2/"$sampleName"_treat_pileup.bdg > $outdir2/"$sampleName".sort.bdg
	awk '{if ($1 != "chrEBV") print $0}' $outdir2/"$sampleName".sort.bdg > $outdir2/"$sampleName".tmp.sort.bdg
	$CUTRUNPATH/bedGraphToBigWig $outdir2/"$sampleName".tmp.sort.bdg $chromsizedir/hg38.chrom.sizes $outdir2/"$sampleName".sorted.bw
	rm -rf $outdir2/"$sampleName".sort.bdg
	rm -rf $outdir2/"$sampleName".tmpsort.bdg

done
echo "Finished"
date