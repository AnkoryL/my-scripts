input <- c("/data/AnkoryL46/results/counts")
library(limma)
library(edgeR)

counts_data <- read.table(paste0(input, "/", "counts_rsem_after_star.txt"))
coldata <- read.table(paste0(input, "/","status.tsv"))
colnames(coldata) <- c("sample", "teg")
m <- coldata[,2]
design_matrix <- model.matrix(~m)
rownames(design_matrix) <- colnames(counts_data)
f <- factor(design_matrix[,2], levels = c("0", "1"))

d0 <- DGEList(counts_data)
d0 <- calcNormFactors(d0)
keep <- rowSums(d0$counts >= 10) >= 3
d <- d0[keep,]

mmatrix <- model.matrix(~0 + f)
y <- voom(d, mmatrix, plot = T)
fit <- lmFit(y, mmatrix)
contr <- makeContrasts(f1 - f0, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)
l <- length(which(top.table$adj.P.Val <= 0.05))
top.table$Gene <- rownames(top.table)
top.table1 <- top.table[1:l,c("Gene", names(top.table)[1:6])]

library("biomaRt")
library("openxlsx")
mart_gene_names = useDataset("mmusculus_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
gene_names_query <- getBM(values = rownames(top.table1), filters = ("ensembl_gene_id"),attributes = c("ensembl_gene_id","external_gene_name","description"), mart = mart_gene_names)
gene_all <- data.frame(gene_names_query, top.table1[2], top.table1[5], top.table1[6])
#write.table(file = "/data/AnkoryL46/results/degs/limmavoom_degs.txt", gene_all, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.xlsx(gene_all, file = "/data/AnkoryL46/results/degs/limmavoom_degs.xlsx", asTable = FALSE)

