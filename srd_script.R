library(ShortRead)

input <- c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw")
output <- c("/data/AnkoryL46/srd_output")

file_list <- list.files(input)
setwd(input)


for (i in 1:length(file_list)) {
  fq <- readFastq(file_list[i])
  reads_median <- median(width(sread(fq)))
  reads_var <- var(width(sread(fq)))
  reads_mean <- mean(width(sread(fq)))
  
}

fq <- readFastq(your_output.truncated)
  reads_median <- median(width(sread(fq)))
  reads_var <- var(width(sread(fq)))
  reads_mean <- mean(width(sread(fq)))


/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw/210618_HSGA.Slominskiy_mRNA_2021.slo19.fastq
/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw/210618_HSGA.Slominskiy_mRNA_2021.slo17.fastq


AdapterRemoval --trimqualities --file1 210618_HSGA.Slominskiy_mRNA_2021.slo19.fastq
AdapterRemoval --trimns --file1 210618_HSGA.Slominskiy_mRNA_2021.slo19.fastq
AdapterRemoval --file1 210618_HSGA.Slominskiy_mRNA_2021.slo19.fastq --basename output_slo19 --trimns --trimqualities

/data/AnkoryL46/a_r/210618_HSGA.Slominskiy_mRNA_2021.slo19.fastq