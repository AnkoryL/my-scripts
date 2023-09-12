cd /data/ivan/tools/STAR/bin/Linux_x86_64
 ./STAR 

#для индексных файлов
cd /data/AnkoryL46/star/input_ref
gunzip *.gz
#system(paste0("gunzip ", input_ref, "/*.gz")))

cd /data/ivan/tools/STAR/bin/Linux_x86_64
./STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir /data/AnkoryL46/star/indices \
--genomeSAindexNbases 12 \
--outFileNamePrefix /data/AnkoryL46/star/indices \
--genomeFastaFiles /data/AnkoryL46/star/input_ref/Mus_musculus.GRCm39.dna.toplevel.fa \
--sjdbGTFfile /data/AnkoryL46/star/input_ref/Mus_musculus.GRCm39.107.gtf

#для картирования
./STAR --genomeDir /data/AnkoryL46/star/indices \
--runThreadN 16 \
--readFilesIn /data/AnkoryL46/a_r_out/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncated \
--outFileNamePrefix /data/AnkoryL46/star/results/210527_HSGA.Slominsky_mRNA_2021.Slo1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within 

#для скриптов
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

 /data/ivan/tools/STAR/bin/Linux_x86_64/STAR \
  --genomeDir /data/AnkoryL46/indices  \
  --outSAMunmapped Within  \
  --outFilterType BySJout  \
  --outSAMattributes NH HI AS NM MD \
  --outFilterMultimapNmax 20  \
  --outFilterMismatchNmax 999  \
  --outFilterMismatchNoverLmax 0.04  \
  --alignIntronMin 20  \
  --alignIntronMax 1000000  \
  --alignMatesGapMax 1000000  \
  --alignSJoverhangMin 8  \
  --alignSJDBoverhangMin 1  \
  --sjdbScore 1  \
  --runThreadN 16  \
  --genomeLoad NoSharedMemory  \
  --outSAMtype BAM Unsorted  \
  --quantMode TranscriptomeSAM  \
  --outSAMheaderHD @HD VN:1.4 SO:unsorted  \
  --outFileNamePrefix /data/AnkoryL46/results/rsem_out_test/Slo1_var2.temp/Slo1_var2  \
  --readFilesIn /data/AnkoryL46/results/trimmed/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncated 

rsem-tbam2gbam /data/AnkoryL46/indices/rsem /data/AnkoryL46/results/rsem_out_test/Slo1_var2..transcript.bam /data/AnkoryL46/results/rsem_out_test/Slo1_var2..genome.bam
rsem-run-gibbs /data/AnkoryL46/indices/rsem /data/AnkoryL46/results/rsem_out_test/Slo1_var2..temp/Slo1_var2. /data/AnkoryL46/results/rsem_out_test/Slo1_var2..stat/Slo1_var2. 200 1000 1 -p 32
rsem-calculate-credibility-intervals /data/AnkoryL46/indices/rsem /data/AnkoryL46/results/rsem_out_test/Slo1_var2..temp/Slo1_var2. /data/AnkoryL46/results/rsem_out_test/Slo1_var2..stat/Slo1_var2. 0.95 1000 50 30000 -p 32

#____________________________________________________________
setwd(star)
for (i in 1:length(trim_file_list)) {
  str <- paste("./STAR --genomeDir ",  index,
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
  " --outFileNamePrefix /data/AnkoryL46/results/rsem_out_test/Slo1_var2.temp/Slo1_var2 ",
  " --readFilesIn /data/AnkoryL46/results/trimmed/210527_HSGA.Slominsky_mRNA_2021.Slo1.fastq.truncated ")
  
system(str)s
}
