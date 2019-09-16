#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -o logs/shuffle_overlap_o
#PBS -e logs/shuffle_overlap_e
#PBS -N shuffle_overlap_te

#1=working directory
#2=annotation file name (eg: exons.bed)
#3=job number (1 to 10)
#4=number of iterations
#5=condition (TB,L12,S12,control)

module load mugqic/bedtools/2.26.0
module load xz/5.2.2
module load gcc/6.1.0
module load intel64/15.3.187
module load tiff/4.0.3
module load bioinformatics/R/3.3.0

cp $1/counts.families.txt $1/counts/counts.families.$3.$2

cp $1/counts.big.families.txt $1/counts_big/counts.big.families.$3.$2 

cp $1/counts.classes.txt $1/counts_classes/counts.classes.$3.$2 

for i in $(seq "$4")
do
        echo $i

#Shuffling peaks
bedtools shuffle -i $1/peaks_by_annotation/peaks.$2 -g $1/genome/hg19.chrom.sizes -incl $1/annotations/$2 > $1/shuffled/shuffled.$3.$2
line=$(head -1 $1/shuffled/shuffled.$3.$2)
echo $line

#Intersecting TEs with shuffled peaks
bedtools intersect -wa -u -a $1/repeats/repeatmasker.mod.bed -b $1/shuffled/shuffled.$3.$2 > $1/overlapping/te.in.shuffled.$3.$2
n=$(wc -l $1/overlapping/te.in.shuffled.$3.$2)
echo $n

echo "Running R scripts"

### Count number of TEs overlapping peaks for each TE family 
Rscript /home/lbogdan/scripts/shuffle/shuffle_overlap_te.R "$1" "$2" "$3" "$i" "families" "counts"

### Count for each big TE family
Rscript /home/lbogdan/scripts/shuffle/shuffle_overlap_te.R "$1" "$2" "$3" "$i" "big.families" "counts_big"

### Count for each TE class
Rscript /home/lbogdan/scripts/shuffle/shuffle_overlap_te.R "$1" "$2" "$3" "$i" "classes" "counts_classes"

done

cp -r $1/counts $HOME/output/$5.results
cp -r $1/counts_big $HOME/output/$5.results
cp -r $1/counts_classes $HOME/output/$5.results

echo "End of script"