suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(biomaRt))
source("/data/ivan/rscript/a_few_common_functions.r")

input_table<-c("/data/AnkoryL46/RNA_seq/mice_mptp_rsem_counts.txt")
count_matrix<-read.table(input_table,sep="\t",quote="",header=TRUE)


input_tags<-c("/data/AnkoryL46/tags.csv")
tag_matrix_new<-read.table(input_tags,sep="\t",header=TRUE)

output_path_local_annotation<-c("/data/ivan/RNA_seq/biomart/mmusculus_biomart_annotation_21_08_2023_unique.csv")
local_ensembl_annotation<-read.table(output_path_local_annotation,sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="")


cn_cm<-colnames(count_matrix)
tm_names<-tag_matrix_new[,1]
rownames(tag_matrix_new)<-tm_names

cn_cm<-gsub("\\.fastq","",cn_cm)
cn_cm<-gsub("\\.1","",cn_cm)
cn_cm<-gsub("\\.truncated","",cn_cm)

any(rownames(tag_matrix_new) != colnames(count_matrix))
# test_df<-data.frame(rownames(tag_matrix_new),colnames(count_matrix))
# test_df

colnames(count_matrix)<-cn_cm
rownames(tag_matrix_new)<-cn_cm

tag_matrix_new<-tag_matrix_new[order(rownames(tag_matrix_new)),]

count_matrix<-count_matrix[,order(colnames(count_matrix))]


any(rownames(tag_matrix_new) != colnames(count_matrix))

# local_ensembl_annotation[,1]
# rownames(count_matrix)
common_gene_names<-intersect(local_ensembl_annotation[,1],rownames(count_matrix))


count_matrix<-count_matrix[rownames(count_matrix) %in% common_gene_names,]
local_ensembl_annotation<-local_ensembl_annotation[local_ensembl_annotation[,1] %in% common_gene_names,]


count_matrix<-count_matrix[order(rownames(count_matrix)),]
local_ensembl_annotation<-local_ensembl_annotation[order(local_ensembl_annotation[,1]),]


any(rownames(count_matrix) != local_ensembl_annotation[,1])

filter_matrix<-apply(count_matrix,1,median)
filter_for_count_matrix<-filter_matrix >= 3

subset_count_matrix_filtered<-count_matrix[filter_for_count_matrix,]
local_ensembl_annotation<-local_ensembl_annotation[filter_for_count_matrix,]

# tag_matrix_new
mptpt<-factor(tag_matrix_new[,3],levels=c("nacl","mptp"))
time<-as.numeric(tag_matrix_new[,4])
set<-factor(tag_matrix_new[,5])
design_matrix<-model.matrix(~mptpt + time)
# design_matrix<-model.matrix(~mptpt)
# design_matrix
rownames(design_matrix)<-colnames(count_matrix)
# design_matrix

output_path<-c("/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_mptp.csv")
output_path_degs<-c("/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_degs_mptp.csv")
out_plot_path<-c("/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_mptp.png")


do_limma(count_matrix=subset_count_matrix_filtered,design_matrix=design_matrix,output_file_result=output_path,
         output_file_degs=output_path_degs,local_ensembl_annotation=local_ensembl_annotation,contrast_num=2,threshold_fc=1.5,threshold_pval=0.05,method_string="fdr",output_plot_path=out_plot_path)

output_patht<-c("/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_mptp_with_time.csv")
output_path_degst<-c("/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_degs_mptp_with_time.csv")

do_limma(count_matrix=subset_count_matrix_filtered,design_matrix=design_matrix,output_file_result=output_patht,
         output_file_degs=output_path_degst,local_ensembl_annotation=local_ensembl_annotation,contrast_num=3,threshold_fc=1.5,threshold_pval=0.05,method_string="fdr",output_plot_path=out_plot_path)
