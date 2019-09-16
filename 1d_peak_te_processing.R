### R SCRIPTS FOR COMPUTING ENRICHMENT IN ATAC-Seq ###

library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(multtest)

setwd("Path/to/dir/Output/S12.output/")

### START of merge_shuffled_counts.R ###

main_dir <- "Path/to/dir/Output/S12.output/shuffle/shuffle_final/counts"
setwd(main_dir)

rep.families <- as.character((read.table("Path/to/dir/Output/te_overlap/rep.families.txt"))[,1])
n.families <- c(1:length(rep.families))

#f.empty.counts.families <- paste0(main_dir, "/counts.families.txt")
df.empty.counts.families <- read.table("Path/to/dir/Output/te_overlap/counts.families.txt")
df.empty.counts.families <- cbind(df.empty.counts.families, as.data.frame(matrix(as.numeric(0), nrow=nrow(df.empty.counts.families), ncol=1)))
names(df.empty.counts.families) <- c("TE.family", "counts.families", "Iterations", "Avg.counts.families")

f.list <- list.files(pattern="bed", full.names=T, recursive=FALSE)
n.files <- length(f.list)

# use regex to remove "./" before file name
f.names <- lapply(f.list, function(x){
  sub(".*./(.*)","\\1",x)
})
f.names <- unlist(f.names)

annotations <- c("exons", "introns", "tss", "promoters", "proximal", "distal", "5UTR", "desert")

f.all.names <- paste0("df.counts.families.all.", annotations, ".txt")

# Create empty count files to fill with the loop below
df.counts.families.all.exons <- df.empty.counts.families
write.table(df.counts.families.all.exons, "df.counts.families.all.exons.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.tss <- df.empty.counts.families
write.table(df.counts.families.all.tss, "df.counts.families.all.tss.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.promoters <- df.empty.counts.families
write.table(df.counts.families.all.promoters, "df.counts.families.all.promoters.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.proximal <- df.empty.counts.families
write.table(df.counts.families.all.proximal, "df.counts.families.all.proximal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.distal <- df.empty.counts.families
write.table(df.counts.families.all.distal, "df.counts.families.all.distal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.5UTR <- df.empty.counts.families
write.table(df.counts.families.all.5UTR, "df.counts.families.all.5UTR.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.desert <- df.empty.counts.families
write.table(df.counts.families.all.desert, "df.counts.families.all.desert.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.families.all.introns <- df.empty.counts.families
write.table(df.counts.families.all.introns , "df.counts.families.all.introns.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

# Loop through counts.families by annotation separately (eg: deal with the 10 files with exon counts.families, then with the 10 files with promoter counts.families, etc.)
for (i in annotations) {
  print(i)
  f.names.sub <- f.names[grep(i, f.names)]
  f.all.name <- f.all.names[grep(i, f.all.names)]
  df.all <- read.table(f.all.name) # table with all counts.families across jobs for one annotation
  # Loop through each TE family
  for (fam in rep.families) {
    counts.families = 0 #contains ALL counts.families for that family across one annotation
    iterations = 0 #should be = 1000 (total number of iterations for one annotation)
    #print(fam)
    # Loop through each of 10 files within the same annotation
    for (file in f.names.sub){
      print(file)
      df <- read.table(file)
      new.counts.families <- df[which(df[,1] == fam), 2]
      new.iterations <- df[which(df[,1] == fam), 3]
      counts.families = counts.families + new.counts.families
      iterations = iterations + new.iterations
    }
    # Add sum of counts.families and iterations in final annotation-specific table 
    df.all[which(df.all[,1] == fam), 2] <- counts.families
    df.all[which(df.all[,1] == fam), 3] <- iterations
    # Take the average of the counts.families for each TE family -> create new column and append # Avg.counts.families = counts.families/Iteration
    df.all[which(df.all[,1] == fam), 4] <- counts.families/iterations
  }
  write.table(df.all, f.all.name, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}

# Combine all counts.families across annotations, one TE family at a time -> use column #4 (avg taken)

df.expected.overlaps <- df.empty.counts.families[,c(1,2)]

for (fam in rep.families) {
  counts.families = 0
  for (file in f.all.names) {
    df <- read.table(file)
    new.counts.families <- df[which(df[,1] == fam), 4]
    counts.families = counts.families + new.counts.families
  }
  df.expected.overlaps[which(df.expected.overlaps[,1] == fam), 2] <- round(counts.families)
}
#Path/to/dir/Output/S12.output/shuffle/shuffle_final/counts
write.table(df.expected.overlaps, "S12_expected_te_overlaps.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### END of merge_shuffled_counts.big.families.R ###


### COMPUTE P VALUES (separating by FAMILY) ###

setwd("Path/to/dir/Output/S12.output")

df.empty.counts.families <- read.table("Path/to/dir/Output/te_overlap/counts.families.txt")

df.rep.families <- (read.table("Path/to/dir/Output/te_overlap/rep.families.txt", stringsAsFactors = F))
rep.families <- df.rep.families[,1]
n.families <- c(1:length(rep.families))

#te.in.peak <- list.files(pattern="te.in.peak.bed", full.names=T, recursive=FALSE)
df.te.in.peak <- read.table("te.in.peak.bed")
#names(df.te.in.peak) <- c("chr","start","end","genoLeft","strand","repName","repClass", "repFamily", "repStart", "repEnd", "repLeft", "id", "bin", "swScore", "milliDiv", "milliDel", "milliIns")

n.te.in.peak <- nrow(df.te.in.peak)
View(df.te.in.peak)

df.observed.overlaps <- df.empty.counts.families[,c(1,2)]

# Create data frame with counts.families by TE family for observed overlaps
for (fam in rep.families) {
  counts.families <- length(which(df.te.in.peak[,5] == fam)) 
  df.observed.overlaps[(which(df.observed.overlaps[,1] == fam)), 2] <- counts.families
}

# Extract data frame with counts.families by TE family for expected (shuffled) overlaps
df.expected.overlaps <- read.table("Path/to/dir/Output/S12.output/shuffle/shuffle_final/counts/S12_expected_te_overlaps.txt")

# Create data frame to place results from binomial test
te.enrichment <- as.data.frame(matrix(as.numeric(0), nrow=length(rep.families), ncol=7))
names(te.enrichment) <- c("TE.family", "Observed.TE.peaks", "Expected.TE.peaks", "Fold.enrichment", "Prob.success", "P.value", "Adj.p.value")
te.enrichment[,1] <- as.data.frame(rep.families)

# Run binomial test
for (fam in rep.families) {
  obs <- df.observed.overlaps[(which(df.observed.overlaps[,1] == fam)), 2]
  exp <- df.expected.overlaps[(which(df.expected.overlaps[,1] == fam)), 2]
  n.te.tot <- df.rep.families[which(df.rep.families[,1]== fam), 2]
  
  # If expected count = 0, set it to 1 to get a fold-enrichment and p-value that make sense
  if (exp==0) {
    exp=1
  }
 
  p <- exp/n.te.tot
  # perform test
  p.val <- as.numeric(binom.test(obs, n.te.tot, p, alternative="greater")$p.value)
  # correct for multiple hypothesis testing
  adj.p.val <- p.val*length(n.families)
  fold <- obs/exp
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
View(te.enrichment)

write.table(te.enrichment, "S12_te_enrichment_families.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

