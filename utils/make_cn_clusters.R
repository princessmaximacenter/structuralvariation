## Cluster copy number profiles with CNpare
## last update 2022-03-09
# Uses genome bins of 1mb like rest of our analyses
# export plots of distance metrix and hclust result 
# incl elbow plot to show # clusters vs  Total (within) sum of squares for a cluster


suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(pheatmap)
  library(RColorBrewer)
  library(CNpare)
  
})



wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.cna.R"))
source(paste0(script_dir,"functions.expression.R"))

## settings
expected_autosomal_length = 3088269832-156040895-57227415


source(paste0(wdir,"wilms.conf"))

results_dir="~/results/wilms_v3/cnpare/"
analysis_type="cnpare_genome_bins."


cohort_sup_table_path="~/data/metadata/cohort.wilms_v3.supplementary_table_1.tsv"

cohort = read.table(cohort_sup_table_path,sep = "\t", header=T) 


#use some templating for this
analysis_type="cnpare_genome_bins.${cancer_type}."

cn_similarity_bins_path_template = paste0(results_dir,analysis_type,"similarities.tsv")
distance_plots_path_template=paste0(results_dir,analysis_type,"distances.pdf")
dendrogram_plots_path_template=paste0(results_dir,analysis_type,"dendrogram.pdf")
cn_similarity_cluster_labels_manhattan_path_template=paste0(results_dir,analysis_type,"cluster_labels_order.manhattan.tsv")
cn_similarity_cluster_labels_euclidean_path_template=paste0(results_dir,analysis_type,"cluster_labels_order.euclidean.tsv")

euclidean_bins_cluster_result_path_template=paste0(results_dir,analysis_type,"clusters_euclidean.tsv")
manhattan_bins_cluster_result_path_template=paste0(results_dir,analysis_type,"clusters_manhattan.tsv")

## annotation ----
heatmap_anno_cols = c("cancer_type","fga")
heatmap_annotation = full_cohort[,heatmap_anno_cols] %>% as.data.frame()
rownames(heatmap_annotation) = full_cohort$patient_label
colnames(heatmap_annotation) = heatmap_anno_cols

## RESOURCES ----

map_template_vars=c('${resources_dir}'=resources_dir)

chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband

chr_arms_df = as.data.frame(chr_arms)


# Genome bins ----

genome_bins = get_genome_windows(target_width=1e6,fixed=F)
genome_bins$bin_id=paste0("bin_",1:length(genome_bins))  
names(genome_bins) = genome_bins$bin_id

genome_bins_metadata= mcols(genome_bins) %>% as.data.frame()

chr_bin_arms_overlaps = get_reciprocal_overlap_pairs(genome_bins,chr_arms,reciprocal_overlap = 0,svtype_matching = F)
bin_chr_arms = chr_bin_arms_overlaps[c("set1","set2")] %>% dplyr::rename(bin_id=set1,chr_arm=set2) %>%
  group_by(bin_id) %>% summarize(chr_arm=paste0(unique(sort(chr_arm)),collapse="_"))

genome_bins_metadata = genome_bins_metadata %>% left_join(bin_chr_arms)
mcols(genome_bins) = genome_bins_metadata
genome_bins$bin_id = paste0(genome_bins$chr_arm,"_",genome_bins$bin_id)
names(genome_bins)=genome_bins$bin_id

bins_to_chr = genome_bins %>% as.data.frame() %>% dplyr::select(seqnames,bin_id,chr_arm) %>% unique()

genome_bins_df=as.data.frame(genome_bins)


# Per cancer type ----
cancer_type_selection="wilms_v3"
cancer_type_selection_vars=c('${cancer_type}'=cancer_type_selection)
  
cn_similarity_bins_path = stri_replace_all_fixed(cn_similarity_bins_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)
distance_plots_path = stri_replace_all_fixed(distance_plots_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)
dendrogram_plots_path = stri_replace_all_fixed(dendrogram_plots_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)

euclidean_bins_cluster_result_path = stri_replace_all_fixed(euclidean_bins_cluster_result_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)
manhattan_bins_cluster_result_path = stri_replace_all_fixed(manhattan_bins_cluster_result_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)


# Read segments ----


cna_seg_file_ext = ".modelFinal.seg"

segments_df = data.frame()

for(pid in cohort$patient_id) {
  patient = cohort %>% filter(patient_id == pid)
 
  run_patient_id=pid
   #load chromosome bands
  #cna_seg_file_ext = ".tumor.modelFinal.seg"
  segments_df_path = paste0(cna_data_dir,patient$tumor_id,"*_WGS*",cna_seg_file_ext)
  segments_df_file=Sys.glob(segments_df_path)
  if(length(segments_df_file)!=1){ 
      print("WARNING: input unclear")
      print(segments_df_path)
      next()
  }
  
  contig_lengths = get_contig_lengths(segments_df_file)
  
  if(sum(contig_lengths$length) != expected_autosomal_length) {
    print("WARNING: autosomal length different than expected, is it really hg38?")
  }
  
  modeled_seg= read_modeled_seg(segments_df_file)
  modeled_seg$from_coordinate = paste0(modeled_seg$seqnames,":",modeled_seg$start,"-",modeled_seg$end)
  
  
  #default
  #modeled_seg_cols=c("cna_id","width","cr_l2fc_50","maf_50","cr_l2fc_10","maf_10","cr_l2fc_90","maf_90")  
  modeled_seg_cols=c(modeled_seg_cols,"from_coordinate")
  modeled_seg_cols = modeled_seg_cols[modeled_seg_cols %in% names(modeled_seg)]
  
  if(exists("run_patient_id")) { modeled_seg$patient_id=run_patient_id }
  modeled_seg$tumor_id=patient$tumor_id
  modeled_seg$patient_label=patient$patient_label
  
  segments_df = rbind(segments_df,modeled_seg)
}



# parse segments

segments_df = segments_df %>% call_cna()
segments_df$patient_cna_id=paste0(segments_df$patient_id,"_",segments_df$cna_id)

segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id


#segments to bins 

bin_cna_overlaps = get_reciprocal_overlap_pairs(segments,genome_bins,reciprocal_overlap = 0,svtype_matching = F)
bin_cna_overlaps = bin_cna_overlaps %>% dplyr::rename( patient_cna_id =set1,bin_id=set2)

segments_df = segments_df %>% left_join(bin_cna_overlaps)

segments_by_bins = get_overlaps_modeled_seg_summary(segments_df,group_cols=c("patient_label","bin_id","seqnames"))

segments_by_bins$bin_id = factor(segments_by_bins$bin_id,levels=unique(genome_bins[order(genome_bins)]$bin_id))

segments_by_bins$seqnames = factor(segments_by_bins$seqnames,
                                   levels=unique(segments_by_bins[gtools::mixedorder(segments_by_bins$seqnames),c("seqnames")]))



#remove gvar/stalk/acen  
segments_by_bins = segments_by_bins %>% filter(!grepl("acen|gvar|stalk",bin_id))

#convert to CN
segments_by_bins = segments_by_bins %>% dplyr::mutate(segVal = (2^(cr_l2fc_50))*2) #reverse log2 and reverse fold change 


## Run CNpare ----
matrix_segments_by_bins = segments_by_bins %>% 
  select(patient_label,bin_id,segVal) %>% unique() %>% pivot_wider(names_from=patient_label,values_from=segVal) %>% as.data.frame()

rownames(matrix_segments_by_bins)=matrix_segments_by_bins$cytoband
matrix_segments_by_bins = as.matrix(matrix_segments_by_bins %>% dplyr::select(-bin_id))

bins_similarities<-getSimilarities(dat1=matrix_segments_by_bins, dat2=matrix_segments_by_bins, method="all",pvalue = T)

cn_similarity_bins = bins_similarities %>% dplyr::rename(patient1 = fileid,patient2=id)

write.table(cn_similarity_bins,cn_similarity_bins_path,quote=F,sep="\t",col.names = T,row.names = F)



## Get clusters semi automated ----
if(FALSE){
  
  
  cn_similarity_bins_path = stri_replace_all_fixed(cn_similarity_bins_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)
  dendrogram_plots_path = stri_replace_all_fixed(dendrogram_plots_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)
  
  euclidean_bins_cluster_result_path = stri_replace_all_fixed(euclidean_bins_cluster_result_path_template,names(cancer_type_selection_vars),cancer_type_selection_vars,vectorize=F)
 
  cn_similarity_bins = read.table(cn_similarity_bins_path,sep="\t",header=T)
  
  euclidean_bins = cn_similarity_bins %>% dplyr::select(patient1,patient2,euclidean) %>% unique() %>%
    pivot_wider(names_from="patient2",values_from="euclidean") %>% as.data.frame()
  rownames(euclidean_bins)=euclidean_bins$patient1
  
  euclidean_bins_matrix=euclidean_bins %>% dplyr::select(-patient1)
  euclidean_bins_matrix = as.matrix(euclidean_bins_matrix)
  euclidean_bins_dist = as.dist(euclidean_bins_matrix)
  
 
  
  euclidean_bins_cluster = hclust(euclidean_bins_dist,method = "ward.D2")

  
  dev.off()
  pdf(dendrogram_plots_path, width=10, height=8)
  plot(euclidean_bins_cluster)
  tree = as.dendrogram(euclidean_bins_cluster)
  
  tree = tree %>% rev()
  plot(tree)
  
  euclidean_bins_cutree = cutree(euclidean_bins_cluster, k = 2:8)
  
  
  wss_clusters_euclidean = apply(euclidean_bins_cutree, 2, function(g) cluster_within_sum_square(euclidean_bins_matrix, g))
  
  wss_clusters_euclidean = as.data.frame(wss_clusters_euclidean) %>% dplyr::rename(wss=wss_clusters_euclidean)
  wss_clusters_euclidean$cluster_cnt = rownames(wss_clusters_euclidean)
  
  p=ggplot(wss_clusters_euclidean,aes(x=cluster_cnt,y=wss)) + geom_line(group="clusters") + ggtitle("Elbow plot for CN similary clusters euclidean") + theme_bw()
  print(p)
  
  
 
  #3 clusters was best
  cluster_k=3
  
  euclidean_bins_cluster_result = euclidean_bins_cutree[,(cluster_k-1)] %>% as.data.frame()
  #euclidean_bins_cluster_result$patient_id=rownames(euclidean_bins_cluster_result)
  euclidean_bins_cluster_result$patient_label=rownames(euclidean_bins_cluster_result)
  euclidean_bins_cluster_result$cluster_id=paste0("CN",euclidean_bins_cluster_result$.)
  euclidean_bins_cluster_result = euclidean_bins_cluster_result[,c("patient_label","cluster_id")] %>% left_join(cohort[,c("patient_label","tumor_id")]) %>% arrange(cluster_id)
  
  write.table(euclidean_bins_cluster_result,euclidean_bins_cluster_result_path,
              quote=F,sep="\t",col.names = T,row.names = F)
  
  
  cluster_k=cn_cluster_manhattan[cancer_type_selection]
  
  
}


