### R SCRIPTS FOR PROCESSING PEAK-TEs ###

library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(multtest)

### START of merge_shuffled_counts.R ###

main_dir <- "Path/to/dir/Output/S12.output/counts_big"
setwd(main_dir)

rep.big.families <- as.character((read.table("Path/to/dir/Output/te_overlap/rep.big.families.txt"))[,1])
n.big.families <- c(1:length(rep.big.families))

df.empty.counts.big.families <- read.table("Path/to/dir/Output/te_overlap/counts.big.families.txt")
df.empty.counts.big.families <- cbind(df.empty.counts.big.families, as.data.frame(matrix(as.numeric(0), nrow=nrow(df.empty.counts.big.families), ncol=1)))
names(df.empty.counts.big.families) <- c("TE.family", "counts.big.families", "Iterations", "Avg.counts.big.families")

f.list <- list.files(pattern="bed", full.names=T, recursive=FALSE)
n.files <- length(f.list)

# use regex to remove "./" before file name
f.names <- lapply(f.list, function(x){
  sub(".*./(.*)","\\1",x)
})
f.names <- unlist(f.names)

annotations <- c("exons", "introns", "tss", "promoters", "proximal", "distal", "5UTR", "desert")

f.all.names <- paste0("df.counts.big.families.all.", annotations, ".txt")

# Create empty count files to fill with the loop below
df.counts.big.families.all.exons <- df.empty.counts.big.families
write.table(df.counts.big.families.all.exons, "df.counts.big.families.all.exons.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.tss <- df.empty.counts.big.families
write.table(df.counts.big.families.all.tss, "df.counts.big.families.all.tss.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.promoters <- df.empty.counts.big.families
write.table(df.counts.big.families.all.promoters, "df.counts.big.families.all.promoters.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.proximal <- df.empty.counts.big.families
write.table(df.counts.big.families.all.proximal, "df.counts.big.families.all.proximal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.distal <- df.empty.counts.big.families
write.table(df.counts.big.families.all.distal, "df.counts.big.families.all.distal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.5UTR <- df.empty.counts.big.families
write.table(df.counts.big.families.all.5UTR, "df.counts.big.families.all.5UTR.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.desert <- df.empty.counts.big.families
write.table(df.counts.big.families.all.desert, "df.counts.big.families.all.desert.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.big.families.all.introns <- df.empty.counts.big.families
write.table(df.counts.big.families.all.introns , "df.counts.big.families.all.introns.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

# Loop through counts.big.families by annotation separately (eg: deal with the 10 files with exon counts.big.families, then with the 10 files with promoter counts.big.families, etc.)
for (i in annotations) {
  print(i)
  f.names.sub <- f.names[grep(i, f.names)]
  f.all.name <- f.all.names[grep(i, f.all.names)]
  df.all <- read.table(f.all.name) # table with all counts.big.families across jobs for one annotation
  # Loop through each TE family
  for (fam in rep.big.families) {
    counts.big.families = 0 #contains ALL counts.big.families for that family across one annotation
    iterations = 0 #should be = 1000 (total number of iterations for one annotation)
    #print(fam)
    # Loop through each of 10 files within the same annotation
    for (file in f.names.sub){
      print(file)
      df <- read.table(file)
      new.counts.big.families <- df[which(df[,1] == fam), 2]
      new.iterations <- df[which(df[,1] == fam), 3]
      counts.big.families = counts.big.families + new.counts.big.families
      iterations = iterations + new.iterations
    }
    # Add sum of counts.big.families and iterations in final annotation-specific table 
    df.all[which(df.all[,1] == fam), 2] <- counts.big.families
    df.all[which(df.all[,1] == fam), 3] <- iterations
    # Take the average of the counts.big.families for each TE family -> create new column and append # Avg.counts.big.families = counts.big.families/Iteration
    df.all[which(df.all[,1] == fam), 4] <- counts.big.families/iterations
  }
  write.table(df.all, f.all.name, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}

# Combine all counts.big.families across annotations, one TE family at a time -> use column #4 (avg taken)

df.expected.overlaps <- df.empty.counts.big.families[,c(1,2)]

for (fam in rep.big.families) {
  counts.big.families = 0
  for (file in f.all.names) {
    df <- read.table(file)
    new.counts.big.families <- df[which(df[,1] == fam), 4]
    counts.big.families = counts.big.families + new.counts.big.families
  }
  df.expected.overlaps[which(df.expected.overlaps[,1] == fam), 2] <- round(counts.big.families)
}
#Path/to/dir/Output/S12.output/shuffle/shuffle_final/counts
write.table(df.expected.overlaps, "S12_expected_te_overlaps.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### END of merge_shuffled_counts.big.big.families.R ###


### COMPUTE P VALUES (separating by FAMILY) ###

setwd("Path/to/dir/Output/S12.output")

df.empty.counts.big.families <- read.table("Path/to/dir/Output/te_overlap/counts.big.families.txt")

df.rep.big.families <- read.table("Path/to/dir/Output/te_overlap/rep.big.families.txt", stringsAsFactors = F)
rep.big.families <- as.character((df.rep.big.families)[,1])
n.big.families <- c(1:length(rep.big.families))


df.te.in.peak <- read.table("te.in.peak.bed")
#names(df.te.in.peak) <- c("chr","start","end","genoLeft","strand","repName","repClass", "repFamily", "repStart", "repEnd", "repLeft", "id", "bin", "swScore", "milliDiv", "milliDel", "milliIns")

n.te.in.peak <- nrow(df.te.in.peak)
View(df.te.in.peak)

df.observed.overlaps <- df.empty.counts.big.families[,c(1,2)]

# Create data frame with counts.big.families by TE family for observed overlaps
for (fam in rep.big.families) {
  counts.big.families <- length(which(df.te.in.peak[,7] == fam)) ############### change this parameter depending on which family level is being evaluated
  df.observed.overlaps[(which(df.observed.overlaps[,1] == fam)), 2] <- counts.big.families
}


# Extract data frame with counts.big.families by TE family for expected (shuffled) overlaps
# Create data frame to place results from binomial test
te.enrichment <- as.data.frame(matrix(as.numeric(0), nrow=length(rep.big.families), ncol=7))
names(te.enrichment) <- c("TE.family", "Observed.TE.peaks", "Expected.TE.peaks", "Fold.enrichment", "Prob.success", "P.value", "Adj.p.value")
te.enrichment[,1] <- as.data.frame(rep.big.families)

# Run binomial test
for (fam in rep.big.families) {
  obs <- df.observed.overlaps[(which(df.observed.overlaps[,1] == fam)), 2]
  exp <- df.expected.overlaps[(which(df.expected.overlaps[,1] == fam)), 2]
  n.te.tot <- df.rep.big.families[which(df.rep.big.families[,1]== fam), 2]
  
  # If expected count = 0, set it to 1 to get a fold-enrichment and p-value that make sense
  if (exp==0) {
    exp=1
  }
 
  p <- exp/n.te.tot
  # perform test
  p.val <- as.numeric(binom.test(obs, n.te.tot, p, alternative="two.sided")$p.value)
  # correct for multiple hypothesis testing
  adj.p.val <- p.val*length(n.big.families)
  #adj.p.val <- as.numeric(p.adjust(p.val,method="bonferroni",n=length(n.big.families)))
  fold <- obs/exp
  #
  new.row <- c(fam, obs, exp, fold, p, p.val, as.numeric(adj.p.val))
  te.enrichment[(which(te.enrichment[,1] == fam)), ] = new.row
} 

te.enrichment$Adj.p.value <- as.numeric(te.enrichment$Adj.p.value)
te.enrichment$Observed.TE.peaks<- as.numeric(te.enrichment$Observed.TE.peaks)
te.enrichment$Expected.TE.peaks<- as.numeric(te.enrichment$Expected.TE.peaks)
p.vals <- te.enrichment$P.value
adj.p.vals <- as.numeric(p.adjust(p.vals,method="BH"))
te.enrichment$Adj.p.value <- adj.p.vals

te.enrichment <- te.enrichment[order(te.enrichment$Adj.p.value),]


write.table(te.enrichment, "S12_te_enrichment_big.families.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

