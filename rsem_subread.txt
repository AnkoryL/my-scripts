setwd("/data/AnkoryL46/indices_rsubread")
buildindex(basename="/data/AnkoryL46/indices_rsubread/Mus_musculus.GRCm39.dna",reference="/data/AnkoryL46/references/Mus_musculus.GRCm39.dna.toplevel.fa")


#----------------------------------------
sub_out <- paste0(output, "/rsubread")
setwd(trim_output)

align(
  
  #  индекс для эталонных последовательностей 
  "/data/AnkoryL46/indices_rsubread/Mus_musculus.GRCm39.dna",
  
  # ввод ридов и вывод 
  readfile1 = trim_file_list[1],
  readfile2 = NULL,
  type = "rna",
  input_format = "FASTQ",
  output_format = "BAM",
  output_file = paste0(sub_out, "/", readfile1, ".subread.", output_format),
  
  # значение смещения, добавляемое к показателям качества Phred баз чтения 
  phredOffset  =  33 ,
  
  # пороги для сопоставления 
  nsubreads  =  10 , 
  TH1  =  3 , 
  TH2  =  1 , 
  maxMismatches  =  3 , 
  
  # уникальное сопоставление и множественное сопоставление 
  unique  =  FALSE , 
  nBestLocations  =  16 , 
  
  # обнаружение indel 
  indels  =  5 , 
  complexIndels  =  FALSE , 
  
  # обрезка чтения 
  nTrim5  =  0 , 
  nTrim3  =  0 , 
  
  # расстояние и ориентация парных концевых считываний
  minFragLength  =  50 , 
  maxFragLength  =  600 , 
  PE_orientation  =  "fr" , 
  
  # количество потоков процессора 
  nthreads  =  16 , 
  
  # группа чтения 
  readGroupID  =  NULL , 
  readGroup  =  NULL , 
  
  # порядок чтения 
  keepReadOrder  =  FALSE , 
  sortReadsByCoordinates  =  FALSE , 
  
  # цветовое пространство читает 
  color2base  =  FALSE , 
  
  # динамическое программирование 
  DP_GapOpenPenalty  =  -1 ,
  DP_GapExtPenalty  =  0 , 
  DP_MismatchPenalty  =  0 , 
  DP_MatchScore  =  2 , 
  
  # определение структурных вариантов 
  detectSV  =  FALSE , 
  
  # аннотация гена 
  useAnnotation  =  TRUE , 
  annot.ext  =  "/data/AnkoryL46/references/Mus_musculus.GRCm39.107.gtf" , 
  isGTF  =  TRUE , 
  GTF.featureType  =  "transcript" , 
  GTF.attrType  =  "gene_id" , 
  chrAliases  =  NULL )

RsubreadUsersGuide()


#--------------------------------------------------------

library(ShortRead)
setwd(trim_output)

system(paste0("rsem-sam-validator /data/AnkoryL46/results/rsubread/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncated.subread.BAM"))
system(paste0("convert-sam-for-rsem /data/AnkoryL46/results/rsubread/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncated.subread.BAM /data/AnkoryL46/results/rsubread/210527_HSGA.Slominsky_mRNA_2021.Slo1.subread.out"))



#----------------------------------------------------------


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
