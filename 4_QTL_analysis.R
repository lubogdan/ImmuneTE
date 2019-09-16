###################################################################
################# Convert reQTL snp IDs to coordinates ############
###################################################################

condition <- "S12"

# Get list of top SNPs corresponding to each DEG
qtl.list <- read.table("Path/to/dir/Output/QTL/nedelec_reqtl.txt", header = T, fill=T, stringsAsFactors = FALSE)

if (condition == "L12") {
  col.sub <- c(1,2,3,5,7,9) 
} else if (condition == "S12") {
  col.sub <- c(1,2,4,6,8,10)
}

qtl.list.sub <- qtl.list[, col.sub]

# remove reQTLs with pval > 0.05
qtl.list.sub.2 <- qtl.list.sub[which(qtl.list.sub[,6] == 1),]

### Convert SNP (rs) ids into genomic coordinates ###

common.snp.150 <- read.table("Path/to/dir/Output/QTL/common.snp150.hg19.bed.txt", fill=T, stringsAsFactors = FALSE)
# remove unwanted chromosomes (eg: hap, chrM)
common.snp.150.sub <- common.snp.150[-c(grep("gl",common.snp.150$V1), grep("hap",common.snp.150$V1), grep("chrM",common.snp.150$V1)),]

qtl.ids <- qtl.list.sub.2[,3]

# For each QTL, get row in SNP coordinates matrix by rsID
snp.coord.rows <- lapply (qtl.ids, function(id){
  print(id)
  which(common.snp.150.sub$V4 == id)
}) 

 # View problematic QTL ids (not present in common snp coordinates)
prob.qtl <- qtl.ids[which(snp.coord.rows == "integer(0)")]

lapply(prob.qtl, function(id){
  qtl.list.sub.2[which(qtl.list.sub.2[,3] == id),]
})

### FOR reQTLs and eQTLs ###
for (id in prob.qtl) {
  print(id)
  qtl.list.sub.2 <- qtl.list.sub.2[-which(qtl.list.sub.2[,3] == id),]
}

# Get SNP coordinates
snp.coord <- common.snp.150.sub[unlist(snp.coord.rows),]

# add genomic coordinates to QTL matrix
qtl.list.sub.3 <- cbind(qtl.list.sub.2, snp.coord)
qtl.list.sub.4 <- qtl.list.sub.3[,c(7:10,1:5)]

# make sure IDs from the two merged matrices are the same -> should get integer(0)
which((qtl.list.sub.4[,4] == qtl.list.sub.4[,7]) == FALSE)
View(qtl.list.sub.4)

# export QTLs with coordinates as BED file
# Headers: c("chr", "start", "end", "rsID", "ensembl_gene_ID", "external_gene_name", "FC_L_top_snp_ID", "FC_L_estimate", "FC_L_pvalue")
write.table(qtl.list.sub.4 , "Path/to/dir/Output/QTL/S12/reqtl.coord.S12.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


### Overlap QTLs with TEs using bedtools on CLUSTER ###
### View results ###

qtl.te <- read.table("Path/to/dir/Output/QTL/100kb.reqtl.near.te.peak.S12.bed", fill=T, stringsAsFactors = FALSE)
#View(qtl.te)

qtl.te.pos <- qtl.te[which(qtl.te[,11] != -1),]
names(qtl.te.pos) <- c("proxy.chr", "proxy.start", "proxy.end", "proxy.SNP.id", "original.SNP.id", "inter.SNP.distance", "gene", "QTL.p.val",
                       "proxy.SNP", "TE.chr", "TE.start","TE.end","TE.strand","TE.family","TE.big.family","TE.class")

View(qtl.te.pos)

write.table(qtl.te.pos , "Path/to/dir/Output/QTL/results/100kb.S12.reqtl.TE.txt", row.names=FALSE, col.names=T, sep="\t", quote=FALSE)


################################################
##### Generate LD intervals based on rAggr #####
################################################

# Input reqtl.coord.S12.bed into rAggr: http://raggr.usc.edu (just copy and paste list of rsIDs in input box)
# Settings: 
# 1- 1000 genomes, phase 3
# 2- populations: CEU and YRI
# 3- max distance: 50 kb
# 4- r2 range -> 0.8 - 1.0 (default)

# name file: Path/to/dir/Output/QTL/S12/raggr.out.S12.union.csv
# http://raggr.usc.edu/Results/Results.aspx?jobId=2116619806

### Format rAggr output ###

raggr.out <- read.table("Path/to/dir/Output/QTL/raggr.out.S12.union.50kb.csv", sep=",", stringsAsFactors = F)
View(raggr.out)

# Remove SNPs with indels in both ref and alternate allele columns
strn <- raggr.out[,11]
n.char <- nchar(strn)
ind.ref <- which(n.char > 1)
#
strn <- raggr.out[,12]
n.char <- nchar(strn)
ind.alt <- which(n.char > 1)
#
index <- unique(c(ind.ref, ind.alt))
#
raggr.out.2 <- raggr.out[-index,]

View(raggr.out.2)

# Extract rsIDs for SNP1
snp1 <- raggr.out.2$V1
ids.snp1 <- c()
for (str in snp1) {
  pieces <- (strsplit(str, ":"))[[1]]
  rsid <- pieces[1]
  ids.snp1 <- c(ids.snp1, rsid)
}

# Extract rsIDs for SNP2
snp2 <- raggr.out.2$V8
ids.snp2 <- c()
for (str in snp2) {
  pieces <- (strsplit(str, ":"))[[1]]
  rsid <- pieces[1]
  ids.snp2 <- c(ids.snp2, rsid)
}

# Add "chr" before chromosome number in columns 2 and 9
chr.snp1 <- paste0("chr",raggr.out.2$V2)
chr.snp2 <- paste0("chr",raggr.out.2$V9)

raggr.out.3 <- cbind(raggr.out.2, ids.snp1, chr.snp1, ids.snp2, chr.snp2)

# Extract important columns
raggr.out.4 <- raggr.out.3[,c(22, 10, 10, 21, 20, 3, 3, 19, 18, 15)]
View(raggr.out.4)

##### Get QTL q-val and gene for each snp from Nedelec reQTL matrix #####

qtl.coord <- read.table("Path/to/dir/Output/QTL/reqtl.coord.S12.bed",stringsAsFactors = FALSE)

# Loop through each proxy SNP and add (empty) gene, ID and pval columns
n.prox <- nrow(raggr.out.4)
gene <- (rep("no", n.prox))
ensemblID <- rep("no", n.prox)
QTLpval <- rep("no", n.prox)
raggr.out.5 <- cbind(raggr.out.4, gene, ensemblID, QTLpval, stringsAsFactors=F)
names(raggr.out.5) <- c("chr.snp2", "start.snp2", "end.snp2", "id.snp2", "chr.snp1", "start.snp1", "end.snp1", "id.snp1", "distance", "pop", "gene", "ENSEMBL.ID", "QTL.p.val")

# Get associated genes from Nedelec QTL matrix and combine info
for (matrix.row in (1:n.prox)) {
  id = as.character(raggr.out.5[matrix.row, 8])
  print(id)
  #get row index
  nedelec.row <- which(qtl.coord$V4 == id)
  #if there is more than one gene associated with a SNP
  if (length(nedelec.row) > 1) {
    print(matrix.row)
    raggr.out.5[matrix.row,11] <- paste((qtl.coord[nedelec.row, 6]), collapse = ",")
    raggr.out.5[matrix.row,12] <- paste((qtl.coord[nedelec.row, 5]), collapse = ",")
    raggr.out.5[matrix.row,13] <- paste((qtl.coord[nedelec.row, 9]), collapse = ",")
  } else {
    #fill matrix
    raggr.out.5[matrix.row,c(11,12,13)] <- qtl.coord[nedelec.row, c(6,5,9)]
  }
}
View(raggr.out.5)

### Add QTLs removed by rAggr ###

# number of QTLs kept
raggr.qtls <- unique(raggr.out.5$id.snp1)
length(raggr.qtls)
# find rsIDs of removed SNPs from Nedelec's matrix
all.qtls <- qtl.coord$V4
length(unique(all.qtls))
# number of missing qtls
length(unique(all.qtls)) - length(raggr.qtls)
# get missing qtls
missing.qtls <- unlist(lapply(all.qtls, function(qtl){
  index <- grep(qtl, raggr.qtls)
  if (length(index) != 1) {
    print(qtl)
  }
}))
# loop through each missing QTL, make its own row and add it to the qtl-proxy list
for (qtl in missing.qtls) {
  nedelec.row <- which(qtl.coord$V4 == qtl)
  # make new row
  snp2.info <- qtl.coord[nedelec.row, c(1:4)]
  snp1.info <- snp2.info
  distance <- 0
  pop <- "none"
  gene <- qtl.coord[nedelec.row, 6]
  ensembl <- qtl.coord[nedelec.row, 5]
  pval <- qtl.coord[nedelec.row, 9]
  # combine and add new row
  new.row <- (cbind(snp2.info, snp1.info, distance, pop, gene, ensembl, pval))
  names(new.row) <- c("chr.snp2", "start.snp2", "end.snp2", "id.snp2", "chr.snp1", "start.snp1", "end.snp1", "id.snp1", "distance", "pop", "gene", "ENSEMBL.ID", "QTL.p.val")
  raggr.out.5 <- rbind(raggr.out.5, new.row)
}
View(raggr.out.5)

write.table(raggr.out.5 , "Path/to/dir/Output/QTL/S12/raggr.S12.formatted.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### Get LD-block intervals ###

raggr.out.5 <- read.table("Path/to/dir/Output/QTL/S12/raggr.S12.formatted.bed", stringsAsFactors = F)
names(raggr.out.5) <- c("chr.snp2", "start.snp2", "end.snp2", "id.snp2", "chr.snp1", "start.snp1", "end.snp1", "id.snp1", "distance", "pop", "gene", "ENSEMBL.ID", "QTL.p.val")
View(raggr.out.5)

ori.snps <- unique(raggr.out.5$id.snp1)
genes <- unique(raggr.out.5$ENSEMBL.ID)

LD.matrix <- as.data.frame(matrix(as.numeric(0), nrow=length(genes), ncol=10))
names(LD.matrix) <- c("Gene.ID", "Gene.name", "Original.SNP", "QTL.p.val", "Chr", "Proxy.SNP.1", "LD.start", "Proxy.SNP.2", "LD.end","LD.width")
LD.matrix[,1] <- genes

for (gene in genes) {
  print(gene)
  #get row to fill
  row.LD <- which(LD.matrix$Gene.ID == gene)
  #select all rows with the gene id
  rows <- which(raggr.out.5$ENSEMBL.ID == gene)
  print(rows)
  #gene.name <- unique(raggr.out.5[rows,7])
  #fill gene name, original snp, pval, chr
  LD.matrix[row.LD,c(2,3,4,5)] <- unique(raggr.out.5[rows,c(11,8,13,1)])
  # get all proxy snps
  #proxy.ids <- raggr.out.5[rows, 4]
  # get all proxy snp starts
  proxy.ids <- raggr.out.5[rows, 2]
  # most upstream snp
  LD.matrix[row.LD,7] <- min(proxy.ids)
  # most downstream snp
  LD.matrix[row.LD,9] <- max(proxy.ids)
  # get block width
  LD.matrix[row.LD,10] <- max(proxy.ids) - min(proxy.ids)
}

LD.blocks <- LD.matrix[,c(5,7,9,1,2,3,4,6,8,10)]
View(LD.blocks)

# Look at LD block distribution
summary(LD.matrix$LD.width)
sd(LD.matrix$LD.width)
#barplot(LD.matrix[order(LD.matrix$QTL.p.val),10], ylim=c(0,100000))
barplot(LD.matrix$LD.width)
# view by decreasing pval
num.pval <- as.numeric(LD.matrix$QTL.p.val)
LD.matrix.2 <- cbind(LD.matrix, num.pval)
LD.matrix.3 <- LD.matrix.2[-which(is.na(num.pval)),c(10,11)]
barplot(LD.matrix.3[order(LD.matrix.3$num.pval), 1], ylim=c(0,100000), ylab="Interval width (bp)", xlab="Increasing QTL pvalue")

# Save small and large LD blocks separately
#large blocks
LD.blocks.large <- LD.blocks[which(LD.blocks$LD.width >= 1000),]
write.table(LD.blocks.large , "Path/to/dir/Output/QTL/S12/S12.LD.blocks.large.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

#small blocks - set coordinates to original reQTLS before saving
LD.blocks.small <- LD.blocks[which(LD.blocks$LD.width < 1000),]
#qtls.small.blocks <- LD.blocks.small$Original.SNP
index <- c(1:nrow(LD.blocks.small))
for (i in index) {
  id <- LD.blocks.small[i,6]
  coord.row <- which(qtl.coord$V4 == id)
  if (length(coord.row > 1)) {
    coord.row <- coord.row[1]
  }
  LD.blocks.small[i,c(1:3)] <- qtl.coord[coord.row, c(1:3)]
}
write.table(LD.blocks.small , "Path/to/dir/Output/QTL/S12/S12.LD.blocks.small.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### enlarge small intervals to 1kb with bedtools slop! ###
# scp Path/to/dir/Output/QTL/S12/S12.LD.blocks.small.bed lbogdan@bourque-mp2.ccs.usherbrooke.ca:/mnt/parallel_scratch_mp2_wipe_on_december_2017/bourque/lbogdan/S12.output/
# cd $scratch/S12.output
# bedtools slop -i S12.LD.blocks.small.bed -g genome/hg19.chrom.sizes -b 500 > S12.LD.blocks.enlarged.bed
# scp lbogdan@bourque-mp2.ccs.usherbrooke.ca:/mnt/parallel_scratch_mp2_wipe_on_december_2017/bourque/lbogdan/S12.output/S12.LD.blocks.enlarged.bed Path/to/dir/Output/QTL/S12/

### Combine small all intervals to a width of at least 1kb ###

LD.blocks.large <- read.table("Path/to/dir/Output/QTL/S12/S12.LD.blocks.large.bed", stringsAsFactors = F)
LD.blocks.enlarged <- read.table("Path/to/dir/Output/QTL/S12/S12.LD.blocks.enlarged.bed", stringsAsFactors = F)

LD.blocks.final <- rbind(LD.blocks.large, LD.blocks.enlarged)
write.table(LD.blocks.final , "Path/to/dir/Output/QTL/S12/S12.LD.blocks.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

# transfer to cluster
# scp Path/to/dir/Output/QTL/S12/S12.LD.blocks.bed lbogdan@bourque-mp2.ccs.usherbrooke.ca:/mnt/parallel_scratch_mp2_wipe_on_december_2017/bourque/lbogdan/S12.output/LD/


################################################
##### Compute LD enrichment in ATAC peaks ######
################################################

### START of merge_shuffled_counts.R ###

main_dir <- "Path/to/dir/Output/S12.output/QTL/counts_LD"
setwd(main_dir)

#f.empty.counts.LD <- paste0(main_dir, "/counts.LD.txt")
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
  #print(fam)
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

# to get observed count: wc -l LD/S12.LD.peaks.bed
# or with TEs: wc -l LD/S12.LD.te.peaks.bed
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

View(df.expected.overlaps)
#Path/to/dir/Output/S12.output/shuffle/shuffle_final/counts
write.table(df.expected.overlaps, "S12_expected_LD_overlaps.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### END of merge_shuffled_counts ###
