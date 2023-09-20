y1 <- t(read.table("/data/AnkoryL46/results/counts/counts_rsem_after_star.txt", 
                 header=T, sep="\t", as.is=TRUE, 
                 check.names=F, comment.char="", 
                 row.names=1 ))

y2 <- t(read.table("/data/AnkoryL46/results/counts/counts_rsem_2.txt", 
                 header=T, sep="\t", as.is=TRUE, 
                 check.names=FALSE, comment.char="", 
                 row.names=1 ))

result_s = cor.test(y1, y2, method = "spearman", exact=FALSE)
result_k = cor.test(y1, y2, method = "kendall", exact=FALSE)

result_s
result_k