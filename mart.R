BiocManager::install("biomaRt")
library("biomaRt")

DEgenes <- read.table('/data/AnkoryL46/results/deseq2_results.txt')

mart_gene_names = useDataset("mmusculus_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))

gene_names_query <- getBM(values = rownames(DEgenes), filters = ("ensembl_gene_id"),attributes = c("ensembl_gene_id","external_gene_name","description"), mart = mart_gene_names)

gene_all <- data.frame(gene_names_query, DEgenes[2], DEgenes[5], DEgenes[6])
  
write.table(file = "/data/AnkoryL46/results/martDEgenes.txt", gene_all, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

