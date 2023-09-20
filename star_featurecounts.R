#folders
maindir <- c("/data/AnkoryL46")
setwd(maindir)

output <- c("/data/AnkoryL46/results")

input <- c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw")
input_ref <- c("/data/AnkoryL46/references")
output_fastqc <- paste0(output, "/fastqc")
trim_output <- paste0(output, "/trimmed")
star_output <- paste0(output, "/star_results")
counts_output <- paste0(output, "/counts")

star <- c("/data/ivan/tools/STAR/bin/Linux_x86_64")
index <- c("/data/AnkoryL46/indices")
counts_file <- c("/data/AnkoryL46/results/counts.txt")

dir.create(output, showWarnings = FALSE)
dir.create(input, showWarnings = FALSE)
dir.create(output_fastqc, showWarnings = FALSE)
dir.create(trim_output, showWarnings = FALSE)
dir.create(star_output, showWarnings = FALSE)
dir.create(counts_output, showWarnings = FALSE)

#fastqc
file_list <- list.files(input)
file_list <- file_list[grep(".fastq", file_list)]
setwd(input)
for (i in 1:length(file_list)) {
  fqc <- paste("fastqc -o ", output_fastqc, file_list[i])
  system(fqc)
}

#adapterremoval
for (i in 1:length(file_list))  {
  trim <- paste("AdapterRemoval --file1 ", file_list[i], " --basename ", trim_output, "/", file_list[i], " --trimns --trimqualities", sep="", collapse=NULL)
  system(trim)
}

#shortread
library(ShortRead)

trim_file_list <- list.files(trim_output)
trim_file_list <- trim_file_list[grep(".fastq.truncated", trim_file_list)]
setwd(trim_output)

for (i in 1:length(trim_file_list)) {
  fq <- readFastq(trim_file_list[i])
  reads_median <- median(width(sread(fq)))
  reads_sd <- sd(width(sread(fq)))
  reads_mean <- mean(width(sread(fq)))
}

#star
setwd(star)
for (i in 1:length(trim_file_list)) {
  str <- paste("./STAR --genomeDir ",  index,
               " --runThreadN 16 ",
               " --readFilesIn ",  trim_output, "/", trim_file_list[i],
               " --outFileNamePrefix ", star_output, "/", trim_file_list[i],
               " --outSAMtype BAM SortedByCoordinate ",
               "--outSAMunmapped Within ", sep="", collapse=NULL)
  system(str)
}

names <- list.files(input, full.names = FALSE)
names <- gsub(".fastq", "", grep(".fastq", names, value = TRUE))

#featurecounts
library("Rsubread")
gtf <- list.files(input_ref, pattern=".*gtf$", full.names=TRUE)
bamFiles <- list.files(star_output, pattern=".*bam$")
setwd(star_output)

for (i in 1:length(bamFiles)) {
  fc <- featureCounts(files = bamFiles[i], GTF.featureType = "exon", GTF.attrType = "gene_id", annot.ext = gtf, nthreads = 16, isGTFAnnotationFile = TRUE)
  if (i == 1) {
    df <- data.frame(fc$counts, stringsAsFactors = FALSE)
  }
  else {
  df <- data.frame(df, fc$counts, stringsAsFactors = FALSE)
  }
}
colnames(df) <- names
write.table(df, file = paste0(counts_output, "/counts.txt"), quote = FALSE, sep = "\t", row.names = T)
