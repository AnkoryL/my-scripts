#folders
maindir <- c("/data/AnkoryL46")
setwd(maindir)

output <- c("/data/AnkoryL46/results")

input <- c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw")
input_ref <- c("/data/AnkoryL46/references")
output_fastqc <- paste0(output, "/fastqc")
trim_output <- paste0(output, "/trimmed")
rsem_output <- paste0(output, "/", "rsem_out_after_star")
star_output <- paste0(output, "/star_out_for_rsem")
counts_output <- paste0(output, "/counts")

star <- c("/data/ivan/tools/STAR/bin/Linux_x86_64")
rsem <- c("/data/ivan/tools/RSEM/rsem-calculate-expression")
index <- c("/data/AnkoryL46/indices")
index_file <- c("/data/AnkoryL46/indices/rsem")

dir.create(output, showWarnings = FALSE)
dir.create(output_fastqc, showWarnings = FALSE)
dir.create(trim_output, showWarnings = FALSE)
dir.create(star_output, showWarnings = FALSE)
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

#star
setwd(star)
trim_file_list <- list.files(trim_output)
trim_file_list <- trim_file_list[grep(".fastq.truncated", trim_file_list)]

for (i in 1:length(trim_file_list)) {
  str <- paste0("./STAR --genomeDir ",  index,
               " --outSAMunmapped Within ",
               " --outFilterType BySJout ",
               " --outSAMattributes NH HI AS NM MD ",
               " --outFilterMultimapNmax 20 ",
               " --outFilterMismatchNmax 999 ",
               " --outFilterMismatchNoverLmax 0.04 ",
               " --alignIntronMin 20 ",
               " --alignIntronMax 1000000 ",
               " --alignMatesGapMax 1000000 ",
               " --alignSJoverhangMin 8 ",
               " --alignSJDBoverhangMin 1 ",
               " --sjdbScore 1 ",
               " --runThreadN 16 ",
               " --genomeLoad NoSharedMemory ",
               " --outSAMtype BAM Unsorted ",
               " --quantMode TranscriptomeSAM ",
               " --outSAMheaderHD @HD VN:1.4 SO:unsorted ",
               " --outFileNamePrefix ", star_output, "/", trim_file_list[i],
               " --readFilesIn ", trim_output, "/", trim_file_list[i])
  system(str)
}

#rsem
library(ShortRead)
setwd(trim_output)

#system(paste0("rsem-sam-validator /data/AnkoryL46/results/star_out_for_rsem/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncatedAligned.toTranscriptome.out.bam"))
#system(paste0("convert-sam-for-rsem /data/AnkoryL46/results/star_out_for_rsem/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncatedAligned.toTranscriptome.out.bam /data/AnkoryL46/results/star_out_for_rsem/210527_HSGA.Slominsky_mRNA_2021.Slo1.out"))

rsem_path <- gsub(".fastq.truncated", "", trim_file_list)
rsem_path <- paste(rsem_output, rsem_path, sep="/", collapse=NULL)
bamFiles <- list.files(star_output, pattern=".*.toTranscriptome.out.bam$", full.names = TRUE)

for (i in 1:length(trim_file_list)) {
  fq <- readFastq(trim_file_list[i])
  reads_sd <- sd(width(sread(fq)))
  reads_mean <- mean(width(sread(fq)))
  f <- paste0(rsem, 
              " --fragment-length-mean ", reads_mean,
              " --fragment-length-sd ", reads_sd,
              " -p 32 --calc-ci --ci-memory 50000 ",
              " --alignments ", bamFiles[i],
              " ", index_file,
              " ", rsem_path[i])
  system(f)
}

#out
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
write.table(input_matrix, file = paste0(counts_output, "/counts_rsem_after_star.txt"), quote = FALSE, sep = "\t", row.names = TRUE)