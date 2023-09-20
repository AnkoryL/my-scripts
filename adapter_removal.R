input <- c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw/")
trim_output <- c("/data/AnkoryL46/a_r_out")

file_list <- list.files(input)
file_list <- file_list[grep(".fastq", file_list)]
setwd(input)

for (i in 1:length(file_list)) {
  fu <- paste("AdapterRemoval --file1 ", file_list[i], " --basename ", trim_output, "/", file_list[i], " --trimns --trimqualities", sep="", collapse=NULL)
  system(fu) 
}
 

library(ShortRead)

tr_input <- c("/data/AnkoryL46/a_r_out/")
read_output <- c("/data/AnkoryL46/srd_output/")

trim_file_list <- list.files(tr_input)
trim_file_list <- trim_file_list[grep(".fastq.truncated", trim_file_list)]
setwd(tr_input)

for (i in 1:length(trim_file_list)) {
  fq <- readFastq(trim_file_list[i])
  reads_median <- median(width(sread(fq)))
  reads_sd <- sd(width(sread(fq)))
  reads_mean <- mean(width(sread(fq)))
  
  print(reads_median)
  print(reads_sd)
  print(reads_mean)
}


