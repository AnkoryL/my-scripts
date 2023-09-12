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

#deseq2
setwd("/data/AnkoryL46/results/counts")
library("DESeq2")

counts_data <- read.table('counts_rsem_after_star.txt')
colData <- read.table('status.tsv')
a <- colData[,2]
design_matrix <- model.matrix(~a)

rownames(design_matrix) <- colnames(counts_data)

dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = design_matrix,
                              design = design_matrix)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)

res_ordered <- res[order(res$pvalue),]
res_sig <- subset(res_ordered, padj <= 0.05 & (log2FoldChange <= -0.58 | log2FoldChange >= 0.58) )
res_sig <- as.data.frame(res_sig)

write.table(res, file = "/data/AnkoryL46/results/degs/deseq2_results_all_genes.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(res_sig, file = "/data/AnkoryL46/results/degs/deseq2_results_degs.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(rownames(res_sig), file = "/data/AnkoryL46/results/degs/deseq2_list_names_degs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

library("biomaRt")
library("openxlsx")
DEgenes <- read.table('/data/AnkoryL46/results/degs/deseq2_results_degs.txt')
mart_gene_names <- read.table("/data/AnkoryL46/mouse_mart_genes.csv",sep="\t",quote="",header=TRUE)
mart_gene_names <- mart_gene_names[(mart_gene_names[,1] %in% rownames(DEgenes)),]
any(mart_gene_names[,1] != rownames(DEgenes))
#DEgenes <- DEgenes[order(rownames(DEgenes)), ]

gene_full <- data.frame(mart_gene_names, DEgenes[2], DEgenes[5], DEgenes[6])
#write.table(file = "/data/AnkoryL46/results/degs/deseq2_degs.txt", gene_full, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.xlsx(gene_full, file = "/data/AnkoryL46/results/degs/deseq2_degs.xlsx", asTable = FALSE)

t <- rowSums(counts(dds) >= 3) >= 10
ddsT <- dds[t,]
write.table(rownames(counts(ddsT)), file = "/data/AnkoryL46/results/degs/deseq2_list_names_genes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#--volcanoplot

mart_gene_names = useDataset("mmusculus_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
gene_names <- getBM(values = rownames(de), filters = ("ensembl_gene_id"),attributes = c("ensembl_gene_id","external_gene_name"), mart = mart_gene_names)
a <- data.frame(row.names(de) == gene_names$ensembl_gene_id)
#------------------------------------
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.0005), col="red")

de$Diffexpressed <- "Not sig"

de$Diffexpressed[de$log2FoldChange > 0.36 & de$padj < 0.0445] <- "Up"
de$Diffexpressed[de$log2FoldChange < -0.5453 & de$padj < 0.048] <- "Down"

de$delabel <- NA

ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=Diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("green", "grey", "red")) +
  geom_vline(xintercept=c(-0.545, 0.35), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#histogram deseq2
genes <- read.xlsx("/data/AnkoryL46/results/degs/deseq2_degs.xlsx")
counts_data <- read.table("/data/AnkoryL46/results/counts/counts_rsem_after_star.txt")
colData <- read.table("/data/AnkoryL46/results/counts/status.tsv")
output <- c("/data/AnkoryL46/results/degs/histograms_deseq2")
output1 <- c("/data/AnkoryL46/results/degs/barplot_deseq2")
dir.create(output)
dir.create(output1)

colData$V1 <- (colData$V1 = str_trunc(colData$V1, 8))
colnames(counts_data) <- colData$V1
counts_data <- counts_data[(rownames(counts_data) %in% genes$ensembl_gene_id),]
any(rownames(counts_data) != genes$ensembl_gene_id)

write.xlsx(counts_data, file = "/data/AnkoryL46/results/degs/count_deseq2_degs.xlsx", rownames = T, colnames = F, asTable = FALSE)

for (i in 1:length(counts_data)) {
  png(paste0(output, "/", genes$external_gene_name[i], ".png"), width = 800, height = 600)
  hist(as.numeric(counts_data[i,]), breaks = 15, xlab = "Counts", main = genes$external_gene_name[i])
  dev.off()
}

for (i in 1:length(counts_data)) {
  png(paste0(output1, "/", genes$external_gene_name[i], ".png"), width = 1200, height = 800)
  barplot(as.numeric(counts_data[i,]), ylab = "Counts", xlab = "Samples", main = genes$external_gene_name[i], names.arg = colnames(counts_data))
  dev.off()
}


