#nohup Rscript /data/AnkoryL46/rscripts/rscript/mice_mptpt_geo/RNA_seq_rsem_mice_mptpt_6m_2m_sets_grcm39.r > /data/AnkoryL46/nohup_logs/RNA_seq_rsem_mice_mptpt_6m_2m_sets_grcm39.r.log 2>&1 &
suppressMessages(library(ShortRead))
source("/data/ivan/rscript/a_few_common_functions.r")

# step 1 - quality control

# load folder

main_output_folder<-c("/data/AnkoryL46/RNA_seq/mice_mptpt_mrna_june25_grcm39")
input_folder<-c("/data/ivan/RNA_seq/mice_mptpt_mrna_june25/raw")
output_folder_trim<-paste0(main_output_folder,"/","fastq_trimmed")
output_folder_qc<-paste0(main_output_folder,"/","rqc")
output_folder_alignments<-paste0(main_output_folder,"/","aligned")

index_file<-c("/data/ivan/references/mouse/grcm39/index/rsem_index")
rsem<-c("/data/ivan/tools/RSEM/rsem-calculate-expression")
star<-c("/data/ivan/tools/STAR-2.7.10b/bin/Linux_x86_64")
output_file_fpkm<-c("/data/AnkoryL46/RNA_seq/mice_mptpt_mrna_june25_grcm39/mice_mptpt_june_25_2022_grcm39_rsem_fpkm.txt")
output_file_counts<-c("/data/AnkoryL46/RNA_seq/mice_mptpt_mrna_june25_grcm39/mice_mptpt_june_25_2022_grcm39_rsem_counts.txt")
# Quality control before truncation. Done with rqc
#rqc(path = input_folder, pattern=".fastq",outdir=output_folder_qc, file="rqc_report_before_trimming", openBrowser = FALSE)



dir.create(main_output_folder, showWarnings = FALSE)
dir.create(input_folder, showWarnings = FALSE)
dir.create(output_folder_qc, showWarnings = FALSE)
dir.create(output_folder_alignments, showWarnings = FALSE)
dir.create(output_folder_trim, showWarnings = FALSE)

filelist<-list.files(input_folder)
filelist<-filelist[grep(".fastq",filelist)]





# trimmed<-TRUE
trimmed<-FALSE


if (!trimmed) {
for (i in 1:length(filelist)) {
com<-paste("AdapterRemoval --file1 ",input_folder,"/",filelist[i]," --basename ",output_folder_trim,"/",filelist[i], " --trimns --trimqualities", sep="",collapse=NULL)
	
# AdapterRemoval --file1 /data/rats_transcriptomes/raw/170322_HSGA.Slominsky_RNA_2017.Slo_G6_II.fastq --basename output_single --trimns --trimqualities
system(com)
}
}
for_alignment<-paste0(output_folder_trim,"/",filelist,".truncated")

done_alignments<-grep(".genes.results",list.files(output_folder_alignments),value=TRUE)
done_alignments<-gsub(".genes.results","",done_alignments)
done_alignments<-gsub(".gz","",done_alignments)
done_alignments<-unique(done_alignments)

alignments<-gsub(".fastq","",filelist)
alignments<-gsub(".truncated","",alignments)
alignments<-gsub(".gz","",alignments)
alignments<-unique(alignments)

alignments_fullpath<-paste(output_folder_alignments,alignments,sep="/",collapse=NULL)

#alignment with RSEM. RSEM uses information about mean length and distribution of reads, so we have to get this data for it
# aligned<-TRUE
aligned<-FALSE
#i<-1
#reads_length_mean<-190.813
#reads_length_sd<-27.70993
if (!aligned) {
	for (i in (1:length(alignments))) {
		if (alignments[i] %in% done_alignments) {next}

		fastq_file<-readFastq(for_alignment[i])
		reads_length_mean<-mean(width(sread(fastq_file)))
		reads_length_sd<-sd(width(sread(fastq_file)))
		fastq_file<-c()
		
#no gzip
		#com<-paste0(rsem," --star --star-path ",star, " \\\n", "--phred33-quals"," \\\n","--fragment-length-mean ", reads_length_mean," \\\n", " --fragment-length-sd ",reads_length_sd ," \\\n", "-p 11 --calc-ci --ci-memory 1024"," \\\n", for_alignment[i]," \\\n",index_file," \\\n" ,alignments_fullpath[i])
#gzip
		#--gzipped-read-file
		com<-paste0(rsem," --output-genome-bam --star --star-path ",
		star, " \\\n ", "--phred33-quals"," \\\n ","--fragment-length-mean ",
		 reads_length_mean," \\\n ", " --fragment-length-sd ",reads_length_sd,
		 " \\\n ", "-p 64 --calc-ci --ci-memory 61024 ", 
		  for_alignment[i]," \\\n "," \\\n ",index_file," \\\n " ,alignments_fullpath[i])
		print(com)
		system(com)
		#rsem-calculate-expression --bowtie-path /sw/bowtie --phred64-quals --fragment-length-mean 150.0 --fragment-length-sd 35.0 -p 8 --output-genome-bam --calc-ci --ci-memory 1024 /data/mmliver.fq /ref/mouse_0 mmliver_single_quals
	}
}
#Now we need to assemble the count matrix for noramlization and processing with limma
filelist_aligned_results<-list.files(output_folder_alignments,full.names=TRUE)
filelist_aligned_results<-grep("genes.result",filelist_aligned_results,value=TRUE)
filelist_aligned_results_names_only<-list.files(output_folder_alignments,full.names=FALSE)
filelist_aligned_results_names_only<-gsub(".genes.results","",grep("genes.result",filelist_aligned_results_names_only,value=TRUE))


input_frame<-data.frame(Gene_id=character(), Transcript_id=character(), Length=numeric())
for (i in 1:length(filelist_aligned_results)) {
	rsem_results<-read.table(filelist_aligned_results[i],sep="\t",quote="",header=TRUE)
		if (i==1) {
		input_matrix<-matrix(,nrow=nrow(rsem_results),ncol=length(filelist_aligned_results))
		rownames(input_matrix)<-rsem_results[,1]
		input_frame<-rsem_results[,1:3]
		}
	input_matrix[,i]<-unlist(rsem_results[,5],use.names=FALSE)
}

colnames(input_matrix)<-filelist_aligned_results_names_only

write.table(input_matrix,output_file_counts,quote = FALSE,sep = "\t",row.names = TRUE)


input_frame<-data.frame(Gene_id=character(), Transcript_id=character(), Length=numeric())
for (i in 1:length(filelist_aligned_results)) {
	rsem_results<-read.table(filelist_aligned_results[i],sep="\t",quote="",header=TRUE)
		if (i==1) {
		input_matrix<-matrix(,nrow=nrow(rsem_results),ncol=length(filelist_aligned_results))
		rownames(input_matrix)<-rsem_results[,1]
		input_frame<-rsem_results[,1:3]
		}
	input_matrix[,i]<-unlist(rsem_results[,7],use.names=FALSE)
}

colnames(input_matrix)<-filelist_aligned_results_names_only

write.table(input_matrix,output_file_fpkm,quote = FALSE,sep = "\t",row.names = TRUE)



