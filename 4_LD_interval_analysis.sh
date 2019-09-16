
### IMPORTANT VARIABLE TO DEFINE ###
condition=S12
peaks=

dir1=$scratch/input_files
dir2=$scratch/S12.output

cd $dir2/LD
bedtools sort -i S12.LD.blocks.bed > tmp
mv tmp S12.LD.blocks.bed

cdsc
cd shuffle

while read l; do
echo $l	
for i in {1..10}
do
echo $i
qsub -F "$dir2 $l $i 100 $condition" -o $dir2/logs/$l.$i.o -e $dir2/logs/$l.$i.e -l walltime=24:00:00 shuffle_overlap_LD.sh
done	
done<$dir2/list_files_no_desert.txt 

# submit separate script to shuffle peaks in gene desert
for i in {1..10}
do
echo $i
qsub -F "$dir2 desert.bed $i 100 $condition" -o $dir2/logs/desert.bed.$i.o -e $dir2/logs/desert.bed.$i.e -l walltime=12:00:00 shuffle_overlap_LD_desert.sh
done

#submit separate script to shuffle peaks in introns
for i in {1..50}
do
echo $i
qsub -F "$dir2 introns.bed $i 20 $condition" -o $dir2/logs/introns.bed.$i.o -e $dir2/logs/introns.bed.$i.e -l walltime=24:00:00 shuffle_overlap_LD.sh
done

#### Bash script ####

#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -o logs/shuffle_overlap_o
#PBS -e logs/shuffle_overlap_e
#PBS -N shuffle_overlap_LD

#1=working directory
#2=annotation file name (eg: exons.bed)
#3=job number (1 to 10)
#4=number of iterations
#5=condition (TB,S12,S12,control)

module load mugqic/bedtools/2.26.0
module load xz/5.2.2
module load gcc/6.1.0
module load intel64/15.3.187
module load tiff/4.0.3
module load bioinformatics/R/3.3.0

cp $1/counts.LD.txt $1/counts_LD/counts.LD.$3.$2

for i in $(seq "$4")
do
        echo $i

#Shuffling peaks
bedtools shuffle -i $1/peaks_by_annotation/peaks.$2 -g $1/genome/hg19.chrom.sizes -incl $1/annotations/$2 > $1/shuffled_LD/shuffled.$3.$2
line=$(head -1 $1/shuffled_LD/shuffled.$3.$2)
echo $line

#Intersecting LD intervals with shuffled peaks
bedtools intersect -wa -u -a $1/shuffled_LD/shuffled.$3.$2 -b $1/LD/$5.LD.blocks.bed > $1/overlapping_LD/LD.in.shuffled.$3.$2
n=$(wc -l $1/overlapping_LD/LD.in.shuffled.$3.$2)
echo $n

echo "Running R scripts"

### Count number of LDs overlapping peaks 
Rscript /home/lbogdan/scripts/shuffle/shuffle_overlap_LD.R "$1" "$2" "$3" "$i" "LD" "counts_LD"
done

echo "End of script"








