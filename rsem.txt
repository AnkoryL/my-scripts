#folders
maindir <- c("/data/AnkoryL46")
setwd(maindir)

output <- c("/data/AnkoryL46/results")

input <- c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw")
input_ref <- c("/data/AnkoryL46/references")
output_fastqc <- paste0(output, "/fastqc")
trim_output <- paste0(output, "/trimmed")
rsem_output <- paste0(output, "/", "rsem_results")
counts_output <- paste0(output, "/counts")

star <- c("/data/ivan/tools/STAR/bin/Linux_x86_64")
rsem <- c("/data/ivan/tools/RSEM/rsem-calculate-expression")
index_file <- c("/data/AnkoryL46/indices/rsem")
counts_file <- c("/data/AnkoryL46/results/counts.txt")

dir.create(output, showWarnings = FALSE)
dir.create(output_fastqc, showWarnings = FALSE)
dir.create(trim_output, showWarnings = FALSE)
dir.create(counts_output, showWarnings = FALSE)
dir.create(rsem_output, showWarnings = FALSE)

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

#rsem
library(ShortRead)
trim_file_list <- list.files(trim_output)
trim_file_list <- trim_file_list[grep(".fastq.truncated", trim_file_list)]
setwd(trim_output)

rsem_path <- gsub(".fastq.truncated", "", trim_file_list)
rsem_path <- paste(rsem_output, rsem_path, sep="/", collapse=NULL)

for (i in 1:length(trim_file_list)) {
    fq <- readFastq(trim_file_list[i])
    reads_sd <- sd(width(sread(fq)))
    reads_mean <- mean(width(sread(fq)))
f <- paste0(rsem, " --output-genome-bam --star --star-path ", star,
                  " --fragment-length-mean ", reads_mean,
                  " --fragment-length-sd ", reads_sd,
                  " -p 32 --calc-ci --ci-memory 30000 ", trim_file_list[i],
                  " ", index_file,
                  " ", alignments_fullpath[i])
system(f)
}

#таблицы

names <- list.files(input, full.names = FALSE)
names <- gsub(".fastq", "", grep(".fastq", names, value = TRUE))
filelist_rsem_output <- list.files(rsem_output, full.names = TRUE)
filelist_rsem_output <- grep("genes.result", filelist_rsem_output, value = TRUE)

for (i in 1:length(filelist_rsem_output)) {
    rsem_out <- read.table(filelist_rsem_output[i], sep = "\t", quote = "", header = TRUE)
  if (i==1) {
    input_matrix <- matrix(, nrow = nrow(rsem_out), ncol = length(filelist_rsem_output))
  }
    input_matrix[,i] <- unlist(rsem_out[,5], use.names = FALSE)
}
rownames(input_matrix) <- rsem_out[,1]
colnames(input_matrix) <- names
write.table(input_matrix, file = paste0(counts_output, "/counts_rsem_2.txt"), quote = FALSE, sep = "\t", row.names = TRUE)
