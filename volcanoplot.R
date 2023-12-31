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
