### R SCRIPTS FOR PROCESSING PEAK-TEs ###

library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(multtest)

### START of merge_shuffled_counts.R ###

main_dir <- "Path/to/dir/Output/L12.output/counts_classes"
setwd(main_dir)

rep.classes <- as.character((read.table("Path/to/dir/Output/te_overlap/rep.classes.txt"))[,1])
n.classes <- c(1:length(rep.classes))

#f.empty.counts.classes <- paste0(main_dir, "/counts.classes.txt")
df.empty.counts.classes <- read.table("Path/to/dir/Output/te_overlap/counts.classes.txt")
df.empty.counts.classes <- cbind(df.empty.counts.classes, as.data.frame(matrix(as.numeric(0), nrow=nrow(df.empty.counts.classes), ncol=1)))
names(df.empty.counts.classes) <- c("TE.family", "counts.classes", "Iterations", "Avg.counts.classes")

f.list <- list.files(pattern="bed", full.names=T, recursive=FALSE)
n.files <- length(f.list)

# use regex to remove "./" before file name
f.names <- lapply(f.list, function(x){
  sub(".*./(.*)","\\1",x)
})
f.names <- unlist(f.names)

annotations <- c("exons", "introns", "tss", "promoters", "proximal", "distal", "5UTR", "desert")

f.all.names <- paste0("df.counts.classes.all.", annotations, ".txt")

# Create empty count files to fill with the loop below
df.counts.classes.all.exons <- df.empty.counts.classes
write.table(df.counts.classes.all.exons, "df.counts.classes.all.exons.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.tss <- df.empty.counts.classes
write.table(df.counts.classes.all.tss, "df.counts.classes.all.tss.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.promoters <- df.empty.counts.classes
write.table(df.counts.classes.all.promoters, "df.counts.classes.all.promoters.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.proximal <- df.empty.counts.classes
write.table(df.counts.classes.all.proximal, "df.counts.classes.all.proximal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.distal <- df.empty.counts.classes
write.table(df.counts.classes.all.distal, "df.counts.classes.all.distal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.5UTR <- df.empty.counts.classes
write.table(df.counts.classes.all.5UTR, "df.counts.classes.all.5UTR.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.desert <- df.empty.counts.classes
write.table(df.counts.classes.all.desert, "df.counts.classes.all.desert.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.classes.all.introns <- df.empty.counts.classes
write.table(df.counts.classes.all.introns , "df.counts.classes.all.introns.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

# Loop through counts.classes by annotation separately (eg: deal with the 10 files with exon counts.classes, then with the 10 files with promoter counts.classes, etc.)
for (i in annotations) {
  print(i)
  f.names.sub <- f.names[grep(i, f.names)]
  f.all.name <- f.all.names[grep(i, f.all.names)]
  df.all <- read.table(f.all.name) # table with all counts.classes across jobs for one annotation
  # Loop through each TE family
  for (fam in rep.classes) {
    counts.classes = 0 #contains ALL counts.classes for that family across one annotation
    iterations = 0 #should be = 1000 (total number of iterations for one annotation)
    #print(fam)
    # Loop through each of 10 files within the same annotation
    for (file in f.names.sub){
      print(file)
      df <- read.table(file)
      new.counts.classes <- df[which(df[,1] == fam), 2]
      new.iterations <- df[which(df[,1] == fam), 3]
      counts.classes = counts.classes + new.counts.classes
      iterations = iterations + new.iterations
    }
    # Add sum of counts.classes and iterations in final annotation-specific table 
    df.all[which(df.all[,1] == fam), 2] <- counts.classes
    df.all[which(df.all[,1] == fam), 3] <- iterations
    # Take the average of the counts.classes for each TE family -> create new column and append # Avg.counts.classes = counts.classes/Iteration
    df.all[which(df.all[,1] == fam), 4] <- counts.classes/iterations
  }
  write.table(df.all, f.all.name, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}

# Combine all counts.classes across annotations, one TE family at a time -> use column #4 (avg taken)

df.expected.overlaps <- df.empty.counts.classes[,c(1,2)]

for (fam in rep.classes) {
  counts.classes = 0
  for (file in f.all.names) {
    df <- read.table(file)
    new.counts.classes <- df[which(df[,1] == fam), 4]
    counts.classes = counts.classes + new.counts.classes
  }
  df.expected.overlaps[which(df.expected.overlaps[,1] == fam), 2] <- round(counts.classes)
}
#Path/to/dir/Output/L12.output/shuffle/shuffle_final/counts
write.table(df.expected.overlaps, "L12_expected_te_overlaps.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### END of merge_shuffled_counts.big.classes.R ###


### COMPUTE P VALUES (separating by FAMILY) ###

setwd("Path/to/dir/Output/L12.output")

df.empty.counts.classes <- read.table("Path/to/dir/Output/te_overlap/counts.classes.txt")

df.rep.classes <- read.table("Path/to/dir/Output/te_overlap/rep.classes.txt", stringsAsFactors = F)
rep.classes <- as.character((df.rep.classes)[,1])
n.classes <- c(1:length(rep.classes))

#te.in.peak <- list.files(pattern="te.in.peak.bed", full.names=T, recursive=FALSE)
df.te.in.peak <- read.table("te.in.peak.bed")
#names(df.te.in.peak) <- c("chr","start","end","genoLeft","strand","repName","repClass", "repFamily", "repStart", "repEnd", "repLeft", "id", "bin", "swScore", "milliDiv", "milliDel", "milliIns")

n.te.in.peak <- nrow(df.te.in.peak)
View(df.te.in.peak)

df.observed.overlaps <- df.empty.counts.classes[,c(1,2)]

# Create data frame with counts.classes by TE family for observed overlaps
for (fam in rep.classes) {
  counts.classes <- length(which(df.te.in.peak[,6] == fam)) ############### change this parameter depending on which family level is being evaluated
  df.observed.overlaps[(which(df.observed.overlaps[,1] == fam)), 2] <- counts.classes
}


# Extract data frame with counts.classes by TE family for expected (shuffled) overlaps

# Create data frame to place results from binomial test
te.enrichment <- as.data.frame(matrix(as.numeric(0), nrow=length(rep.classes), ncol=7))
names(te.enrichment) <- c("TE.family", "Observed.TE.peaks", "Expected.TE.peaks", "Fold.enrichment", "Prob.success", "P.value", "Adj.p.value")
te.enrichment[,1] <- as.data.frame(rep.classes)

# Run binomial test
for (fam in rep.classes) {
  obs <- df.observed.overlaps[(which(df.observed.overlaps[,1] == fam)), 2]
  exp <- df.expected.overlaps[(which(df.expected.overlaps[,1] == fam)), 2]
  n.te.tot <- df.rep.classes[which(df.rep.classes[,1]== fam), 2]
  
  # If expected count = 0, set it to 1 to get a fold-enrichment and p-value that make sense
  if (exp==0) {
    exp=1
  }
 
  p <- exp/n.te.tot
  # perform test
  p.val <- as.numeric(binom.test(obs, n.te.tot, p, alternative="two.sided")$p.value)
  # correct for multiple hypothesis testing
  adj.p.val <- p.val*length(n.classes)
  #adj.p.val <- as.numeric(p.adjust(p.val,method="bonferroni",n=length(n.classes)))
  fold <- obs/exp
  #
  new.row <- c(fam, obs, exp, fold, p, p.val, as.numeric(adj.p.val))
  te.enrichment[(which(te.enrichment[,1] == fam)), ] = new.row
} 

# Organize final table
te.enrichment$Adj.p.value <- as.numeric(te.enrichment$Adj.p.value)
te.enrichment$Observed.TE.peaks<- as.numeric(te.enrichment$Observed.TE.peaks)
te.enrichment$Expected.TE.peaks<- as.numeric(te.enrichment$Expected.TE.peaks)
p.vals <- te.enrichment$P.value
adj.p.vals <- as.numeric(p.adjust(p.vals,method="BH"))
te.enrichment$Adj.p.value <- adj.p.vals
te.enrichment <- te.enrichment[order(te.enrichment$Adj.p.value),]

# Save table
write.table(te.enrichment, "Path/to/dir/Output/L12.output/counts_classes/L12_te_enrichment_classes.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
