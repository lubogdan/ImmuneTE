### Count LD-peak overlaps ###

args <- commandArgs(trailingOnly = TRUE)
main_dir <- as.character(args[1])
anno_file_name <- as.character(args[2])
n.job <- as.character(args[3])
n.i <- as.character(args[4])
fam.type <- as.character(args[5])
out_dir <- as.character(args[6])

setwd(main_dir)

### Define variables, files and paths ###

# Define file names
name.LD.in.shuffled <- paste("LD.in.shuffled", n.job, anno_file_name, sep=".")
name.counts <- paste("counts", fam.type, n.job, anno_file_name, sep=".")

# Get counts file
f.counts <- paste0(main_dir, "/", out_dir, "/", name.counts)
df.counts <- read.table(f.counts)
names(df.counts) <- c(fam.type, "counts", "iterations")
print(df.counts[1,])

# Get file w/ LDs overlapping peaks
print(paste0(main_dir,"/overlapping/",name.LD.in.shuffled))
f.LD.in.shuffled <- paste0(main_dir,"/overlapping/",name.LD.in.shuffled)
# check file size and count number of overlaps
if((file.info(f.LD.in.shuffled))$size == 0){
  new.count = 0
} else {
  df.LD.in.shuffled <- read.table(f.LD.in.shuffled)
  new.count <- nrow(df.LD.in.shuffled)
}

### Add new count to number from previous iterations ###
count = as.numeric(df.counts[1, 2])
count = count + new.count
df.counts[1, 2] = count
df.counts[1, 3] = n.i

print(df.counts[1,])

print("Updated counts table")

# Export new counts file 
write.table(df.counts, f.counts, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

print("Exported new counts table")


