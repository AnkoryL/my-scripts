#FIRST INPUT STARTS HERE

base_path<-"/data/AnkoryL46/RNA_seq/mice_mptp_three_sets/wgcna"
setwd(base_path)

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set

all_file_path<-"/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_with_time.csv"
all_genes_file<-read.table(all_file_path,sep="\t",header=TRUE,quote="")



deg_file_path<-"/data/AnkoryL46/RNA_seq/output_limma/mice_mptp_three_datasets_degs_with_time.csv"
deg_file<-read.table(deg_file_path,sep="\t",header=TRUE,quote="")


# Need to combine 2 read tables into one

input_fpkm_1 = read.table("/data/AnkoryL46/RNA_seq/mice_mptpt_mrna_june25_grcm39/mice_mptpt_june_25_2022_grcm39_rsem_fpkm.txt")
input_fpkm_2 = read.table("/data/AnkoryL46/RNA_seq/mice_mptpt_21_nov_2022/mice_mptpt_21_nov_2022_rsem_fpkm.txt")

# Take a quick look at what is in the data set:

input_raw_1 = read.table("/data/AnkoryL46/RNA_seq/mice_mptpt_mrna_june25_grcm39/mice_mptpt_june_25_2022_grcm39_rsem_counts.txt")
input_raw_2 = read.table("/data/AnkoryL46/RNA_seq/mice_mptpt_21_nov_2022/mice_mptpt_21_nov_2022_rsem_counts.txt")


inputdata_fpkm<-cbind(input_fpkm_1,input_fpkm_2)

inputdata_raw<-cbind(input_raw_1,input_raw_2)

inputdata_merged_fpkm_path = c("/data/AnkoryL46/RNA_seq/mice_mptp_three_sets/wgcna/mice_mptp_rsem_fpkm_raw_merged.txt")
inputdata_merged_raw_path = c("/data/AnkoryL46/RNA_seq/mice_mptp_three_sets/wgcna/mice_mptp_rsem_counts_raw_merged.txt")

# 
write.table(file=inputdata_merged_fpkm_path,sep="\t",x=inputdata_merged_raw_path,col.names=TRUE,row.names=TRUE)

write.table(file=inputdata_merged_raw_path,sep="\t",x=inputdata_raw,col.names=TRUE,row.names=TRUE)


med<-apply(inputdata_raw,1,median)



filter_for_count_matrix<-med>=3

inputdata_fpkm<-inputdata_fpkm[filter_for_count_matrix,]
inputdata_raw<-inputdata_raw[filter_for_count_matrix,]


datExpr0 = as.data.frame(t(inputdata_fpkm));



gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
multiExpr = vector(mode = "list")

datExpr0<-datExpr0[gsg$goodSamples,gsg$goodGenes]
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.905,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#SECOND INPUT STARTS HERE

#PUT SELECTED SOFT THRESHOLD HERE
power_chosen<-26

net = blockwiseModules(datExpr0, power = power_chosen,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "hcmptom",
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];



genelist<-vector(mode="list",length=length(table(moduleLabels))-1)

for (i in 1:(length(table(moduleLabels))-1)) {
  genes_in_module<-names(moduleLabels)[moduleLabels==i]
  genelist[[i]]<-genes_in_module
}

genelist_deg<-vector(mode="list",length=length(genelist))

for (i in 1:length(genelist)) {
  genelist_deg[[i]]<-intersect(genelist[[i]],deg_file[,1])
  
}

newwd<-paste0(base_path,"/","sublists_v_1")

unlink(newwd,recursive = TRUE)
dir.create(newwd, showWarnings = FALSE)
setwd(newwd)

outname<-paste0("module_",1:length(genelist_deg),".csv")
outname_deg<-paste0("module_deg_",1:length(genelist_deg),".csv")

for (i in 1:length(genelist_deg)) {
  out_deg<-deg_file[deg_file[,1] %in% genelist_deg[[i]],]
  out_all<-all_genes_file[all_genes_file[,1] %in% genelist[[i]],]
  
  write.table(x=out_all,file=outname[i],sep="\t",quote=FALSE)
  if (nrow(out_deg)>=1) {
    write.table(x=out_deg,file=outname_deg[i],sep="\t",quote=FALSE)
  }
}



#SECOND INPUT UNTILL HERE
