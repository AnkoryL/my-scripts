
library(WGCNA)
library(org.Mm.eg.db)
library(clusterProfiler)

base_path<-"/data/AnkoryL46/RNA_seq/mice_mptp_three_sets/wgcna"
out_folder_name<-"cluster_profiler_enrichment"
out_folder_path<-paste0(base_path,"/",out_folder_name)
background_path<-"/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_with_time.csv"

# Load the WGCNA package
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set

background_list<-read.table(background_path,sep="\t",header=TRUE,quote="")
background_list<-background_list[,1]

newwd<-paste0(base_path,"/","sublists_v_1")

list_files<-list.files(newwd)
list_paths<-list.files(newwd,full.name=TRUE)

unlink(out_folder_path,recursive = TRUE)
dir.create(out_folder_path, showWarnings = FALSE)

i<-1

outnames<-gsub("\\.csv$","_go_bp_enrichment.csv",list_files)
outpaths<-paste0(out_folder_path,"/",outnames)

for (i in 1:length(list_files)) {
  input_cluster<-read.table(list_paths[i],sep="\t",quote="",header=TRUE)
  if (nrow(input_cluster)<=2) {next}
  tryCatch({
    enrich<-enrichGO(gene=input_cluster[,1],OrgDb=org.Mm.eg.db,keyType = "ENSEMBL",ont = "BP",pAdjustMethod = "fdr",universe=background_list,qvalueCutoff = 0.05)
    out<-enrich@result
    out<-out[out[,9]>=2,]
    if (nrow(out)>=1) {
      out<-out[out[,6] <= 0.05,]
      if (nrow(out)>=1) {
        write.table(file=outpaths[i],x=out,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
      }
    }
  }, warning=function(war) {
    print(paste("Warning caught: ",war))
  }, error=function(err) {
    print(paste("Error caught: ",err))
  })
}
