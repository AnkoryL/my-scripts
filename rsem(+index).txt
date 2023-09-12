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
rsem <- c("/data/ivan/tools/RSEM/rsem-calculate-expression")
index <- c("/data/AnkoryL46/indices/rsem")
rsem_index <- c("/data/ivan/tools/RSEM/rsem-prepare-reference")

setwd("/data/AnkoryL46/indices/")
#index
f <- paste0(rsem_index, " --gtf /data/AnkoryL46/references/Mus_musculus.GRCm39.107.gtf ", 
            "-p 16 ", "--star ", 
            "--star-path ", star,
            " /data/AnkoryL46/references/Mus_musculus.GRCm39.dna.toplevel.fa ",
            "/data/AnkoryL46/indices/rsem")
system(f)

#rsem
  for (i in 1:length(trim_file_list)) {
      fq <- readFastq(trim_file_list[i])
      reads_sd <- sd(width(sread(fq)))
      reads_mean <- mean(width(sread(fq)))      
	f <- paste0(rsem, " --output-genome-bam --star --star-path ", star,
                  " --fragment-length-mean ", reads_mean,
                  " --fragment-length-sd ", reads_sd,
                  " -p 32 --calc-ci --ci-memory 50000 ", trim_file_list[i],
                  " ", index,
                  " ", alignments_fullpath[i])
  system(f)
}

#rsem without star

#1 - star index

f <- paste0(rsem_index, " --gtf /data/AnkoryL46/references/Mus_musculus.GRCm39.107.gtf ", 
            "-p 16 ", "--star ", 
            "--star-path ", star,
            " /data/AnkoryL46/references/Mus_musculus.GRCm39.dna.toplevel.fa ",
            "/data/AnkoryL46/indices/rsem")
system(f)

cd /data/ivan/tools/STAR/bin/Linux_x86_64
./STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir /data/AnkoryL46/star/indices \
--genomeSAindexNbases 12 \
--outFileNamePrefix /data/AnkoryL46/star/indices \
--genomeFastaFiles /data/AnkoryL46/star/input_ref/Mus_musculus.GRCm39.dna.toplevel.fa \
--sjdbGTFfile /data/AnkoryL46/star/input_ref/Mus_musculus.GRCm39.107.gtf

#2 - star again with new index
 
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

#3 - rsem
  for (i in 1:length(trim_file_list)) {
      fq <- readFastq(trim_file_list[i])
      reads_sd <- sd(width(sread(fq)))
      reads_mean <- mean(width(sread(fq)))      
	f <- paste0(rsem, " --output-genome-bam --star --star-path ", star,
                  " --fragment-length-mean ", reads_mean,
                  " --fragment-length-sd ", reads_sd,
                  " -p 32 --calc-ci --ci-memory 50000 ", trim_file_list[i],
                  " ", index,
                  " ", alignments_fullpath[i])
  system(f)
}

