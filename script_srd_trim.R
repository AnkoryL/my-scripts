library(ShortRead)

input <- c("/data/AnkoryL46/a_r_out/")
read_output <- c("/data/AnkoryL46/srd_output/")

file_list <- list.files(input)
file_list <- file_list[grep(".fastq.truncated", file_list)]
setwd(input)

for (i in 1:length(file_list)) {
  fq <- readFastq(file_list[i])
  
  reads_median <- median(width(sread(fq)))
  reads_sd <- sd(width(sread(fq)))
  reads_mean <- mean(width(sread(fq)))
  
  
  print(reads_median)
  print(reads_sd)
  print(reads_mean)
}



