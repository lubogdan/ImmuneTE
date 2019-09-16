
### Compute enrichment of TE families near differentially expressed genes following infection ###

library(multtest)
library(GenomicRanges)
library(ggplot2)


setwd("Path/to/dir/Output/S12.output")

df.deg.control <- read.table("Path/to/dir/Output/distance_to_gene/control_deg.txt")

df.deg.pacis <- read.table("Path/to/dir/Output/distance_to_gene/Pacis_DEG.txt")
names(df.deg.pacis) <- c("Ensembl.id", "Gene.name", "log2FC", "P.value", "FDR")
# keep only genes with FDR < 0.05
df.deg.pacis <- df.deg.pacis[which(df.deg.pacis$FDR < 0.05),]

df.deg.nedelec <- read.table("Path/to/dir/Output/distance_to_gene/Nedelec_DEG.txt")
names(df.deg.nedelec) <- c("Ensembl.id", "Gene.name", "INF.L.logFC", "INF_S_logFC", "INF.L.pvalue", "INF.S.pvalue", "INF.L.fdr", "INF.S.fdr")

# GET DEGs
df.deg.L <- df.deg.nedelec[c(which(df.deg.nedelec$INF.L.logFC < -2), which(df.deg.nedelec$INF.L.logFC > 2)), ]
df.deg.S <- df.deg.nedelec[c(which(df.deg.nedelec$INF_S_logFC < -2), which(df.deg.nedelec$INF_S_logFC  > 2)), ]

df.deg <- df.non.deg.L

### CONVERT ENSEMBL ID TO GENOMIC COORDINATES ###

df.ens.genes <- read.table("Path/to/dir/Output/distance_to_gene/ensembl_genes.bed")
View(df.ens.genes)

deg.ens.ids <- as.character(df.deg$Ensembl.id)

num.transcripts <- nrow(df.ens.genes)

df.deg.coord <- as.data.frame(matrix(as.numeric(0), nrow=num.transcripts, ncol=6))
names(df.deg.coord) <- c("chr", "start", "end", "strand", "gene.id", "transcript.id")

i = 1
for (id in deg.ens.ids) {
  #print(id)
  rows <- which(df.ens.genes$V13==id)
  #print(rows)
  # get cols V3,V5,V6, V13, V2 from df.ens.genes
  for (row in rows) {
    df.deg.coord[i,1] <- as.character(df.ens.genes[row, 3])
    df.deg.coord[i,2] <- (df.ens.genes[row, 5])
    df.deg.coord[i,3] <- (df.ens.genes[row, 6])
    df.deg.coord[i,4] <- as.character(df.ens.genes[row, 4])
    df.deg.coord[i,5] <- as.character(df.ens.genes[row, 13])
    df.deg.coord[i,6] <- as.character(df.ens.genes[row, 2])
    i = i + 1
  }
}

df.deg.coord.raw <- df.deg.coord
df.deg.coord <- df.deg.coord[which(df.deg.coord$chr != "0"),]
write.table(df.deg.coord , "Path/to/dir/Output/distance_to_gene/L.non.deg.coordinates.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


### Add gene names to DEG coordinate data frame ###

df.deg <- read.table("Path/to/dir/Output/distance_to_gene/S.deg.coordinates.txt", stringsAsFactors = F)

df.deg.nedelec <- read.table("Path/to/dir/Output/distance_to_gene/Nedelec_DEG.txt", stringsAsFactors = F)
names(df.deg.nedelec) <- c("Ensembl.id", "Gene.name", "INF.L.logFC", "INF_S_logFC", "INF.L.pvalue", "INF.S.pvalue", "INF.L.fdr", "INF.S.fdr")

rows <- c(1:nrow(df.deg))
gene.names <- c()

for (row in rows){
  id <- df.deg[row, 5]
  gene <- df.deg.nedelec[which(df.deg.nedelec$Ensembl.id == id), 2]
  gene.names <- c(gene.names, gene)
}

df.deg.2 <- cbind(df.deg, gene.names)

#write.table(df.deg.2, "Path/to/dir/Output/distance_to_gene/L.deg.coordinates.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)




### GET FILES ###

rep.families <- as.character((read.table("Path/to/dir/Output/distance_to_gene/rep.families.txt"))[,1])
n.families <- length(rep.families)

# peak-TEs
df.peak.te.all <- read.table("Path/to/dir/Output/S12.output/te.in.peak.bed", stringsAsFactors = F)
df.peak.te.all <- df.peak.te.all[-c(grep("gl",df.peak.te.all$V1), grep("hap",df.peak.te.all$V1), grep("chrM",df.peak.te.all$V1), grep("chrY",df.peak.te.all$V1)),]
gr.peak.te.all <- makeGRangesFromDataFrame(df.peak.te.all, keep.extra.columns=TRUE, seqnames.field="V1", start.field="V2", end.field="V3", strand.field="V4")

# non-peak-TEs
df.non.peak.te.all <- read.table("Path/to/dir/Output/S12.output/te.not.in.peak.bed")
df.non.peak.te.all <- df.non.peak.te.all[-c(grep("gl",df.non.peak.te.all$V1), grep("hap",df.peak.te.all$V1), grep("chrM",df.peak.te.all$V1)),]

# DEGs
df.deg <- read.table("Path/to/dir/Output/distance_to_gene/S.deg.coordinates.txt", stringsAsFactors = F)
gr.deg <- makeGRangesFromDataFrame(df.deg, keep.extra.columns=TRUE, seqnames.field="V1", start.field="V2", end.field="V2", strand.field="V4")

# All genes 
df.ens.genes <- read.table("Path/to/dir/Output/distance_to_gene/ensembl_genes.bed", stringsAsFactors = F)
df.ens.genes <- df.ens.genes[-c(grep("gl",df.ens.genes$V3), grep("hap",df.ens.genes$V3), grep("chrM",df.ens.genes$V3), grep("chrY",df.ens.genes$V3)),]
gr.ens.genes <- makeGRangesFromDataFrame(df.ens.genes, keep.extra.columns=TRUE, seqnames.field="V3", start.field="V5", end.field="V5", strand.field="V4")

### COMPUTE ENRICHMENT ###

# Get TE families that have more than 10 instances in peaks
n.fam = 0
interesting.rep.families <- list()

for (fam in rep.families) {
 
  df.peak.te <- df.peak.te.all[which(df.peak.te.all[,5] == fam),]
  
  if (nrow(df.peak.te) < 10) {
    print(paste("Not enough counts for",fam))
    next
  } else {
    print(paste(nrow(df.peak.te),fam))
    interesting.rep.families <- c(interesting.rep.families,fam)
    n.fam = n.fam + 1
  }
}

interesting.rep.families <- unlist(interesting.rep.families)
n.interesting.families <- length(interesting.rep.families)

# Create table to place p-values
enrichment.near.deg <- as.data.frame(matrix(as.numeric(0), nrow=length(interesting.rep.families), ncol=7))
names(enrichment.near.deg) <- c("TE.family", "Observed.TE.near.DEG", "Expected.TE.near.DEG", "Fold.enrichment", "Prob.success", "P.value", "Adj.p.value")
enrichment.near.deg[,1] <- as.data.frame(interesting.rep.families)

### COMPUTE P-VALUE FOR ENRICHMENT NEAR DEGs FOR EACH TE FAMILIY ###

for (fam in interesting.rep.families) {
  
  df.peak.te <- df.peak.te.all[which(df.peak.te.all[,5] == fam),]
  
  if (nrow(df.peak.te) < 1) {
    print(paste("No counts for",fam))
    next
  } else {
    print(paste(nrow(df.peak.te),fam))
  }
  gr.peak.te <- makeGRangesFromDataFrame(df.peak.te, keep.extra.columns=TRUE, seqnames.field="V1", start.field="V2", end.field="V3", strand.field="V4")
  
  df.non.peak.te <- df.non.peak.te.all[which(df.non.peak.te.all[,5] == fam),]
  gr.non.peak.te <- makeGRangesFromDataFrame(df.non.peak.te, keep.extra.columns=TRUE, seqnames.field="V1", start.field="V2", end.field="V3", strand.field="V4")
  
  ######## OBSERVED DISTRIBUTION ########
  
  hits <- distanceToNearest(gr.peak.te, gr.deg, ignore.strand=T)
  d <- as.data.frame(mcols(hits)$distance)
  
  obs.count <- length(which(d < 100000))
  
  ######## EXPECTED DISTRIBUTION ########
  
  n.peak.te <- nrow(df.peak.te)
  n.non.peak.te <- nrow(df.non.peak.te)
  
  n = 0
  
  all.exp.counts.near.deg = 0
  
  for (i in c(1:1000)) {
    
    # get random sample on non-peak-TEs - sample size = same as peak-TEs
    random <- sample(seq(1,n.non.peak.te,1), n.peak.te)
    
    gr.random <- gr.non.peak.te[random]
    
    hits <- distanceToNearest(gr.random, gr.deg, ignore.strand=T)
    d <- mcols(hits)$distance
    
    count <- length(which(d < 100000))
    all.exp.counts.near.deg <- all.exp.counts.near.deg + count
    
    n = n + 1
    #print(n)
  }
  
  ### COMPUTE P-value #### 
  obs <- obs.count
  exp <- all.exp.counts.near.deg/n
  fold <- obs/exp
  p <- exp/n.peak.te
  # perform test
  p.val <- as.numeric(binom.test(obs, n.peak.te, p, alternative="greater")$p.value)
  # correct for multiple hypothesis testing
  adj.p.val <- as.numeric(p.adjust(p.val,method="BH",n=(n.interesting.families))) 
  #adj.p.val <- as.numeric(p.val*n.interesting.families)
  new.row <- c(fam, obs, exp, fold, p, p.val, as.numeric(adj.p.val))
  enrichment.near.deg[(which(enrichment.near.deg[,1] == fam)), ] = new.row
}

enrichment.near.deg.raw <- enrichment.near.deg

p.vals <- enrichment.near.deg$P.value
adj.p.vals <- as.numeric(p.adjust(p.vals,method="BH"))
enrichment.near.deg$Adj.p.value <- adj.p.vals

enrichment.near.deg <- enrichment.near.deg[order(enrichment.near.deg$Adj.p.value),]

write.table(enrichment.near.deg, "Path/to/dir/Output/S12.output/S12.deg.enrich.gene.100kb.final.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)