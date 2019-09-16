
################################################
##### Compute LD enrichment in ATAC peaks ######
################################################

### START of merge_shuffled_counts.R ###

main_dir <- "Path/to/dir/Output/L12.output/QTL/counts_LD"
setwd(main_dir)

df.empty.counts.LD <- read.table("Path/to/dir/Output/QTL/S12.output/counts.LD.txt")
df.empty.counts.LD <- cbind(df.empty.counts.LD, as.data.frame(matrix(as.numeric(0), nrow=nrow(df.empty.counts.LD), ncol=1)))
names(df.empty.counts.LD) <- c("Interval", "counts.LD", "Iterations", "Avg.counts.LD")

f.list <- list.files(pattern="bed", full.names=T, recursive=FALSE)
n.files <- length(f.list)

# use regex to remove "./" before file name
f.names <- lapply(f.list, function(x){
  sub(".*./(.*)","\\1",x)
})
f.names <- unlist(f.names)

annotations <- c("exons", "introns", "tss", "promoters", "proximal", "distal", "5UTR", "desert")

f.all.names <- paste0("df.counts.LD.all.", annotations, ".txt")

# Create empty count files to fill with the loop below
df.counts.LD.all.exons <- df.empty.counts.LD
write.table(df.counts.LD.all.exons, "df.counts.LD.all.exons.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.tss <- df.empty.counts.LD
write.table(df.counts.LD.all.tss, "df.counts.LD.all.tss.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.promoters <- df.empty.counts.LD
write.table(df.counts.LD.all.promoters, "df.counts.LD.all.promoters.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.proximal <- df.empty.counts.LD
write.table(df.counts.LD.all.proximal, "df.counts.LD.all.proximal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.distal <- df.empty.counts.LD
write.table(df.counts.LD.all.distal, "df.counts.LD.all.distal.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.5UTR <- df.empty.counts.LD
write.table(df.counts.LD.all.5UTR, "df.counts.LD.all.5UTR.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.desert <- df.empty.counts.LD
write.table(df.counts.LD.all.desert, "df.counts.LD.all.desert.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
df.counts.LD.all.introns <- df.empty.counts.LD
write.table(df.counts.LD.all.introns , "df.counts.LD.all.introns.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

# Loop through counts.LD by annotation separately (eg: deal with the 10 files with exon counts.LD, then with the 10 files with promoter counts.LD, etc.)
for (i in annotations) {
  print(i)
  f.names.sub <- f.names[grep(i, f.names)]
  f.all.name <- f.all.names[grep(i, f.all.names)]
  df.all <- read.table(f.all.name) # table with all counts.LD across jobs for one annotation
  #
  counts.LD = 0 #contains ALL counts.LD for that family across one annotation
  iterations = 0 #should be = 1000 (total number of iterations for one annotation)
  # Loop through each of 10 files within the same annotation
  for (file in f.names.sub){
    print(file)
    df <- read.table(file)
    new.counts.LD <- df[, 2]
    new.iterations <- df[, 3]
    counts.LD = counts.LD + new.counts.LD
    iterations = iterations + new.iterations
  }
  # Add sum of counts.LD and iterations in final annotation-specific table 
  df.all[, 2] <- counts.LD
  df.all[, 3] <- iterations
  # Take the average of the counts.LD for each TE family -> create new column and append # Avg.counts.LD = counts.LD/Iteration
  df.all[, 4] <- counts.LD/iterations
  # Export table
  write.table(df.all, f.all.name, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}

# Combine all counts.LD across annotations -> use column #4 (avg taken)

df.expected.overlaps <- df.empty.counts.LD[,c(1,2)]
counts.LD = 0
for (file in f.all.names) {
  df <- read.table(file)
  new.counts.LD <- df[, 4]
  if(is.na(new.counts.LD)){
    new.counts.LD <- 0
  }
  print(new.counts.LD)
  counts.LD = counts.LD + new.counts.LD
}
df.expected.overlaps[, 2] <- round(counts.LD)

### Compute p-value ###

# to get observed count: wc -l LD/L12.LD.peaks.bed
# or with TEs: wc -l LD/L12.LD.te.peaks.bed
obs
exp <- df.expected.overlaps[, 2]
exp
# number of trials = number of TE-peaks
n.te.tot

p <- exp/n.te.tot
# perform test
p.val <- as.numeric(binom.test(obs, n.te.tot, p, alternative="greater")$p.value)
p.val
fold <- obs/exp
fold
#

write.table(df.expected.overlaps, "L12_expected_LD_overlaps.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


