################################################################################
### Full script to run on called peaks for each of the following conditions: ###
################################################################################


### IMPORTANT VARIABLE TO DEFINE ###
condition=S12
peaks=$condition.summits.bed

# Define folders for input and output files
dir1=$HOME/input_files
dir2=$HOME/output/$condition.output

module load xz/5.2.2
module load gcc/6.1.0
module load intel64/15.3.187
module load tiff/4.0.3
module load bioinformatics/R/3.3.0

# Before starting:
### Put narrowPeak file in output folder and call it [condition].peaks.bed ### cp [file.narrowPeak] $dir2/
# cd Path/to/dir/Output/$condition.output

### Create annotation files for shuffling (only need to do this once) - BASH ### Files saved in annotations/
### Modify repeatmasker (only need to do this once) - R ###

# get peak-TEs
cd $dir2
cp $dir1/repeatmasker.mod.bed $dir2/
#remove quotes from file
sed 's/\"//g' resized.$peaks > tmp.bed
mv tmp.bed resized.$peaks
#
bedtools intersect -wa -u -a repeatmasker.mod.bed -b resized.$peaks > te.in.peak.bed
wc -l te.in.peak.bed
# Get TEs NOT overlapping (non-resized) peaks... 
peaks=$condition.peaks.bed
#remove quotes from file
sed 's/\"//g' resized.$peaks > tmp.bed
mv tmp.bed resized.$peaks
bedtools intersect -wa -v -a repeatmasker.mod.bed -b $peaks > te.not.in.peak.bed
# change peaks variable back to summits
peaks=$condition.summits.bed

### Shuffle peaks and intersect with repeatmasker on cluster to get null distribution of counts - BASH and R ###

# On the cluster, copy necessary files to working directory:

cd $dir2
mkdir -p genome
mkdir -p peaks_by_annotation
mkdir -p annotations
mkdir -p repeats
mkdir -p shuffled
mkdir -p overlapping
mkdir -p counts
mkdir -p counts_big
mkdir -p counts_classes

cp $dir1/list_files.txt $dir2/

cp $dir1/list_files_no_desert.txt $dir2/ 

cp $dir1/hg19.chrom.sizes $dir2/genome/

#cp /gs/project/mugqic/projects/lbogdan/annotation_refseq/dedup/final/* $dir2/annotations/
cp $dir1/annotations/* $dir2/annotations/

mv $dir2/repeatmasker.mod.bed $dir2/repeats/repeatmasker.mod.bed

cp $dir1/counts.families.txt $dir2/
cp $dir1/counts.big.families.txt $dir2/
cp $dir1/counts.classes.txt $dir2/

cp $dir1/rep.families.txt $dir2/
cp $dir1/rep.big.families.txt $dir2/
cp $dir1/rep.classes.txt $dir2/

### Separate peaks by annotation ### 

cd $dir2/annotations

for f in *.bed
do
bedtools intersect -u -a ../resized.$peaks -b $f > ../peaks_by_annotation/peaks.$f
done

# get peaks in desert
rm ../peaks_by_annotation/peaks.non.desert.bed
bedtools intersect -v -a ../resized.$peaks -b non.desert.bed > ../peaks_by_annotation/peaks.desert.bed

cd ../peaks_by_annotation/

#desert - (distal + proximal + introns + promoters + tss + exons + 5UTR)
cat peaks.distal.bed peaks.proximal.bed peaks.introns.bed peaks.promoters.bed peaks.tss.bed peaks.exons.bed peaks.5UTR.bed | sort | uniq > tmp.prev.bed
bedtools subtract -a peaks.desert.bed -b tmp.prev.bed > tmp && mv tmp peaks.desert.bed

# distal - (proximal + introns + promoters + tss + exons + 5UTR)
cat peaks.proximal.bed peaks.introns.bed peaks.promoters.bed peaks.tss.bed peaks.exons.bed peaks.5UTR.bed | sort | uniq > tmp.prev.bed
bedtools subtract -a peaks.distal.bed -b tmp.prev.bed > tmp && mv tmp peaks.distal.bed

# proximal - (introns + promoters + tss + exons + 5UTR)
cat peaks.introns.bed peaks.promoters.bed peaks.tss.bed peaks.exons.bed peaks.5UTR.bed | sort | uniq > tmp.prev.bed
bedtools subtract -a peaks.proximal.bed -b tmp.prev.bed > tmp && mv tmp peaks.proximal.bed

# introns - (promoters + tss + exons + 5UTR)
cat peaks.promoters.bed peaks.tss.bed peaks.exons.bed peaks.5UTR.bed | sort | uniq > tmp.prev.bed
bedtools subtract -a peaks.introns.bed -b tmp.prev.bed > tmp && mv tmp peaks.introns.bed

# promoters - (tss + exons + 5UTR)
cat peaks.tss.bed peaks.exons.bed peaks.5UTR.bed | sort | uniq > tmp.prev.bed
bedtools subtract -a peaks.promoters.bed -b tmp.prev.bed > tmp && mv tmp peaks.promoters.bed

# tss - (exons + 5UTR)
cat peaks.exons.bed peaks.5UTR.bed | sort | uniq > tmp.prev.bed
bedtools subtract -a peaks.tss.bed -b tmp.prev.bed > tmp && mv tmp peaks.tss.bed

# exons - 5UTR
bedtools subtract -a peaks.exons.bed -b peaks.5UTR.bed > tmp && mv tmp peaks.exons.bed

rm tmp.prev.bed

# Check number of peaks (should be same as in initial peak file)
n_peaks=0
for f in *.bed
do
echo $f
num=$(wc -l $f | cut -d' ' -f1)
echo $num
n_peaks=$((n_peaks+num))
done
echo $n_peaks
wc -l ../resized.$peaks


### Send jobs with shuffle_overlap_te.sh ###

cdsc
cd shuffle

while read l; do
echo $l	
for i in {1..10}
do
echo $i
qsub -F "$dir2 $l $i 100" -o logs/$condition/$l.$i.o -e logs/$condition/$l.$i.e -l walltime=24:00:00 shuffle_overlap_te.sh
done	
done<$dir2/list_files_no_desert.txt 

# submit separate script to shuffle peaks in gene desert
for i in {1..10}
do
echo $i
qsub -F "$dir2 desert.bed $i" -o logs/$condition/desert.bed.$i.o -e logs/$condition/desert.bed.$i.e -l walltime=12:00:00 shuffle_overlap_te_desert.sh
done

#submit separate script to shuffle peaks in introns
for i in {1..50}
do
echo $i
qsub -F "$dir2 introns.bed $i 20" -o logs/$condition/introns.bed.$i.o -e logs/$condition/introns.bed.$i.e -l walltime=24:00:00 shuffle_overlap_te.sh
done











