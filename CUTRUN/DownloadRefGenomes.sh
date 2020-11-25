#need sample.bed in the current directory

#url is defined as following
url=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFaMasked.tar.gz
#url=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
#http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFaMasked.tar.gz
#http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFaMasked.tar.gz

module load BEDTools

cd assemblies/chrom.hg38
wget $url
tar -zxf hg38.chromFaMasked.tar.gz
ls -1 maskedChroms/chr*.fa|xargs cat > hg38.fa
#should create *.fai index file
bedtools getfasta -fi hg38.fa -bed ../../sample.bed
../../fetchChromSizes hg38 > hg38.chrom.sizes
rm hg38.chromFaMasked.tar.gz
#rm -rf chr*.fa.masked
cd $cur
