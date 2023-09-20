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

res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

res_ordered <- res[order(res$pvalue),]
res_sig <- subset(res_ordered, padj < 0.05)
res_sig <- as.data.frame(res_sig)

write.table(res, file = "/data/AnkoryL46/results/degs/deseq2_results_all_genes.txt", quote = FALSE, sep = "\t", row.names = TRUE)
write.table(res_sig, file = "/data/AnkoryL46/results/degs/deseq2_results_degs.txt", quote = FALSE, sep = "\t", row.names = TRUE)

genelist <- res_sig$padj
names(genelist) <- rownames(res_sig)
write.table(rownames(res_sig), file = "/data/AnkoryL46/results/degs/deseq2_list_names_degs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

library("biomaRt")
library("openxlsx")
DEgenes <- read.table('/data/AnkoryL46/results/degs/deseq2_results_degs.txt')
mart_gene_names <- read.table("/data/AnkoryL46/mouse_mart_genes.csv",sep="\t",quote="",header=TRUE)
mart_gene_names <- mart_gene_names[(mart_gene_names[,1] %in% rownames(DEgenes)),]
any(mart_gene_names[,1] != rownames(DEgenes))
#DEgenes <- DEgenes[order(rownames(DEgenes)), ]

gene_all <- data.frame(mart_gene_names, DEgenes[2], DEgenes[5], DEgenes[6])
#write.table(file = "/data/AnkoryL46/results/degs/deseq2_degs.txt", gene_all, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.xlsx(gene_all, file = "/data/AnkoryL46/results/degs/deseq2_degs.xlsx", asTable = FALSE)

t <- rowSums(counts(dds) >= 3) >= 10
ddsT <- dds[t,]
write.table(rownames(counts(ddsT)), file = "/data/AnkoryL46/results/degs/deseq2_list_names_genes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
