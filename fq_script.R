input <- c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw")
output <- c("/data/AnkoryL46/fq_output")


file_list <- list.files(input)
file_list <- file_list[grep(".fastq", file_list)]
setwd(input)


for (i in 1:length(file_list)) {
  fun_fq <- paste("fastqc -o ", output, file_list[i] )
  
  system(fun_fq) 
}
