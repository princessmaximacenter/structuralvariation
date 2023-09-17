# Annotate merged SVs with repeats
## Reannotate 100 bp distance

suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(dplyr, quietly=TRUE)
  library(GenomicRanges)
  library(openssl)
library(igraph)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
})



root_dir="~/PycharmProjects/"
source(paste0(root_dir,"structuralvariation/sv_functional_analysis/default.conf"))

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.expression.R"))
source(paste0(script_dir,"functions.cna.R"))

source(paste0(wdir,"solids.conf"))

autosomes= c(paste("chr",1:22,sep=""),"chrX")


map_template_vars=c('${resources_dir}'=resources_dir,
                    '${merged_svs_dir}'=merged_svs_dir,
                    #'${cohort_wdir}'=cohort_wdir,
                    #'${cohort_identifier}'=cohort_identifier,
                    '${cna_data_dir}'=cna_data_dir,
                    '${cna_seg_file_ext}'=cna_seg_file_ext,
                    '${sequencing_strategy}'="WGS")


max_repeat_sv_distance=100 #to link svs to repeats

# input
repeatmasker_path = stri_replace_all_fixed(repeatmasker_path_template,names(map_template_vars), map_template_vars,vectorize=F)
segmental_duplications_path = stri_replace_all_fixed(segmental_duplications_path_template,names(map_template_vars), map_template_vars,vectorize=F)

# output 
merged_svs_repeat_path = paste0(cohort_results_dir,"merged_svs.100bp_repeat.20230906.tsv")


# Cohort ----

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")


# Resources ----

## chromosomes ----

chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband

chr_arm_order = chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")]# %>% flatten_chr()
chrom_order = chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")] %>% unique()# %>% flatten_chr()


# Read data ----


## Load SV data ----

unfiltered_svs_df = load_cohort_svs(cohort,c('${merged_svs_dir}'=merged_svs_dir),svs_path_template=svs_union_anno_somatic_multitool_path_template,flag_add_patient_label = T)
unfiltered_svs_df = unfiltered_svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
unfiltered_svs_df = unfiltered_svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
unfiltered_svs_df = annotate_sv_multitool_support(unfiltered_svs_df,sv_tumor_af_threshold = 0.1)
unfiltered_svs_df = unfiltered_svs_df %>% unique()


merged_svs = make_merged_svs(unfiltered_svs_df)
merged_svs = merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  
#export merged svs? > in omics 
#}

merged_svs = merged_svs %>% filter(flag_in_filtered_svs) %>% filter(patient_label %in% cohort$patient_label) %>% unique()
#sets globals therefore no output
resolve_orphan_svs(merged_svs,unfiltered_svs_df)

#check consistency merged svs and svs df 
merged_svs %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% nrow() == 0
unfiltered_svs_df %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% filter(patient_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() == 0

svs_gr = GRanges(svs_df$coordinate)
svs_gr$svtype=svs_df$svtype
#mcols(svs_gr)=svs_df
names(svs_gr)=svs_df$patient_sv_name
  


## Load repeats ----
repeatmasker = read.table(repeatmasker_path,header=T,sep="\t",comment.char = "")
repeatmasker = repeatmasker %>% dplyr::rename("seqnames"="genoName","start"="genoStart","end"="genoEnd")
repeatmasker_bk=repeatmasker
repeatmasker_bk$repClass %>% unique()
repeatmasker = repeatmasker_bk %>% filter(repClass %in% c("LINE","SINE","LTR","Simple_repeat","Low_complexity","Retroposon")) # %>% filter(abs(repLeft)<50)
repeatmasker$repeat_id = paste0("rep_",1:nrow(repeatmasker))
repeatmasker = GRanges(repeatmasker)
names(repeatmasker) = repeatmasker$repeat_id

repeatmasker_df = as.data.frame(repeatmasker)
rownames(repeatmasker_df) = repeatmasker_df$repeat_id


segmental_duplications = read.table(segmental_duplications_path,header=T,sep="\t",comment.char = "")
segmental_duplications = segmental_duplications %>% dplyr::rename("seqnames"="chrom","start"="chromStart","end"="chromEnd","sd_name"="name")
segmental_duplications$segdup_id = paste0("segdup_",1:nrow(segmental_duplications))
segmental_duplications = GRanges(segmental_duplications)
names(segmental_duplications) = segmental_duplications$segdup_id

segmental_duplications_df = as.data.frame(segmental_duplications)


# SVs in repeats ----


#do not resize the svs because we look at breakpoints!!
sv_bp_repeatmasker_overlaps = get_reciprocal_overlap_pairs_start_end(svs_gr,resize_gr_distance(repeatmasker,max_repeat_sv_distance),reciprocal_overlap = 0,svtype_matching = FALSE)
sv_bp_repeatmasker_overlaps = sv_bp_repeatmasker_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, repeat_id=set2) 

sv_bp_repeatmasker_overlaps = sv_bp_repeatmasker_overlaps %>%  
  left_join(repeatmasker_df[,repeatmasker_df_cols] %>% uniq,by=c("repeat_id")) %>%
  left_join(svs_df[,c("patient_sv_name","patient_sv_merged","partner_sv_name")] %>% unique(),by="patient_sv_name")
  
#sv_bp_repeatmasker_overlaps %>% filter(patient_sv_name %in% filter(svs_df,sv_merged_coordinate=="chr13:91046950-91047021:*")$patient_sv_name)

  #use partnering instead of a loop with start/end orientation 
  ## this allows you to annotate CTX on partner just like the the genes
  ## flag if any overlap with repeat
  ## flag if both start/end overlap same repeat family
  
#not filtered on repLeft ...
 #

  sv_bp_repeatmasker_partnered = make_sv_bp_gene_partnered(sv_bp_repeatmasker_overlaps,gene_overlap_cols = c("repClass"),
                                                           patient_sv_col="patient_sv_name",partner_sv_col="partner_sv_name",
                                                           svs_df_cols =c("patient_sv_name","partner_sv_name","svtype"))
  
  sv_bp_repeatmasker_wide = sv_bp_repeatmasker_partnered %>% 
    group_by(patient_sv_name) %>% summarize(start_repeat_class = toString(unique(sort(repClass))),
                                            end_repeat_class = toString(unique(sort(partner_repClass)))) %>% as.data.frame() %>% unique()
  
  sv_bp_repeatmasker_partnered_filtered_repleft = make_sv_bp_gene_partnered(sv_bp_repeatmasker_overlaps %>% filter(abs(repLeft)<50),gene_overlap_cols = c("repClass"),
                                                           patient_sv_col="patient_sv_name",partner_sv_col="partner_sv_name",
                                                           svs_df_cols =c("patient_sv_name","partner_sv_name","svtype"))
  
  sv_bp_repeatmasker_wide_filtered_repleft = sv_bp_repeatmasker_partnered_filtered_repleft %>% 
    group_by(patient_sv_name) %>% summarize(start_repeat_class_filtered_repleft = toString(unique(sort(repClass))),
                                            end_repeat_class_filtered_repleft = toString(unique(sort(partner_repClass)))) %>% as.data.frame() %>% unique()
 
  
  sv_bp_segmental_duplications_overlaps= get_reciprocal_overlap_pairs_start_end(svs_gr,resize_gr_distance(segmental_duplications,max_repeat_sv_distance),reciprocal_overlap = 0,svtype_matching = FALSE)
  
  sv_bp_segmental_duplications_overlaps = sv_bp_segmental_duplications_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, segdup_id=set2) 
  
  sv_bp_segmental_duplications_overlaps = sv_bp_segmental_duplications_overlaps %>% 
      left_join(segmental_duplications_df[,segmental_duplications_df_cols],by=c("segdup_id"))  %>%
    left_join(svs_df[,c("patient_sv_name","patient_sv_merged","partner_sv_name")])#,by="patient_sv_name")
  
  sv_bp_segmental_duplications_partnered = make_sv_bp_gene_partnered(sv_bp_segmental_duplications_overlaps,gene_overlap_cols = c("sd_name"),
                                                                     
                                                           patient_sv_col="patient_sv_name",partner_sv_col="partner_sv_name",
                                                           svs_df_cols =c("patient_sv_name","partner_sv_name","svtype"))
  
## remove annotation if present
unanno_svs_df = svs_df %>% select(-contains("repeat"),-contains("segdup"))
unanno_svs_df = unanno_svs_df %>% left_join(sv_bp_repeatmasker_wide) %>% left_join(sv_bp_repeatmasker_wide_filtered_repleft)
  
  
  col_name_base="flag_repeat"
  unanno_svs_df=flag_sv_overlap(unanno_svs_df,sv_bp_repeatmasker_overlaps,
                           sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                           overlaps_colname=col_name_base)

  
  #both breakpoints same repeat class 
  unanno_svs_df=flag_sv_overlap(unanno_svs_df,
                         sv_bp_repeatmasker_partnered %>% filter(repClass==partner_repClass),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=paste0(c(col_name_base,"both_bp"),collapse = "_"))
  
  col_name_base="flag_repeat_filtered_repleft"
  unanno_svs_df=flag_sv_overlap(unanno_svs_df,sv_bp_repeatmasker_overlaps %>% filter(abs(repLeft)<50),
                           sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                           overlaps_colname=col_name_base)

  
  #both breakpoints same repeat class 
  unanno_svs_df=flag_sv_overlap(unanno_svs_df,
                         sv_bp_repeatmasker_partnered_filtered_repleft %>% filter(repClass==partner_repClass),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=paste0(c(col_name_base,"both_bp"),collapse = "_"))
  
  
  
  col_name_base="flag_segdup"
  unanno_svs_df=flag_sv_overlap(unanno_svs_df,sv_bp_segmental_duplications_overlaps,
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=col_name_base)
  
  unanno_svs_df=flag_sv_overlap(unanno_svs_df,
                         sv_bp_segmental_duplications_partnered %>% filter(!is.na(sd_name) & !is.na(partner_sd_name)),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=paste0(c(col_name_base,"both_bp"),collapse = "_"))
  


## merged svs ----
get_merged_svs_repeat_anno = function(unfiltered_svs_df) {
  merged_svs_anno = unfiltered_svs_df %>%
    group_by(patient_label,svtype,patient_sv_merged,sv_merged_coordinate,flag_in_filtered_svs) %>% summarize(
      flag_repeat_any=any(flag_repeat),flag_repeat_all=all(flag_repeat),
      flag_repeat_both_bp_any=any(flag_repeat_both_bp),flag_repeat_both_bp_all=all(flag_repeat_both_bp),
      flag_segdup_any=any(flag_segdup),flag_segdup_all=all(flag_segdup),
      flag_segdup_both_bp_any=any(flag_segdup_both_bp),flag_segdup_both_bp_all=all(flag_segdup_both_bp),
      start_repeat_class=lst_str(start_repeat_class),end_repeat_class = lst_str(end_repeat_class),
      
      flag_repeat_any_filtered_repleft=any(flag_repeat_filtered_repleft),flag_repeat_all_filtered_repleft=all(flag_repeat_filtered_repleft),
      flag_repeat_both_bp_filtered_repleft_any=any(flag_repeat_filtered_repleft_both_bp),flag_repeat_both_bp_filtered_repleft_all=all(flag_repeat_filtered_repleft_both_bp),
      start_repeat_class_filtered_repleft=lst_str(start_repeat_class_filtered_repleft),end_repeat_class_filtered_repleft = lst_str(end_repeat_class_filtered_repleft)   
      )

   
  return(merged_svs_anno)
}

merged_svs_anno = get_merged_svs_repeat_anno(unanno_svs_df)

merged_svs_anno[duplicated(merged_svs_anno$patient_sv_merged),]

# Export ----

write.table(merged_svs_anno,merged_svs_repeat_path,sep="\t",row.names = F,col.names = T)


