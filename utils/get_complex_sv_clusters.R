## Complex SV clusters 
## Last update: 2023-09


suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  #library(stringdist, quietly=TRUE)
})


wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(wdir,"functions.cna.R"))
source(paste0(wdir,"functions.general.R"))
source(paste0(wdir,"functions.svs.R"))

source(paste0(wdir,"solids.conf"))

# Settings ----
run_date="20230701"
override_clusters=F

amplicons_footprint_size_min = 50e3 # min size of amplicon to be included in analysis
amplicons_sv_distance_max = 1e3 #to link svs to amplicons
complex_sv_bp_distance_max = 5e6

autosomes= c(paste("chr",1:22,sep=""),"chrX")

merged_svs_cols = c("patient_sv_merged","sv_merged_coordinate","svtype","tumor_af_mean","svlen_mean",
                    "cytoband","partner_cytoband","chrom","partner_chrom","chr_arm","partner_chr_arm")


## Cohort paths ----


complex_svs_results_dir=paste0(cohort_results_dir,"complex_sv_clusters/") 


## Inputs
manual_labels_path=paste0(complex_svs_results_dir,"complex_sv.manual_labels.txt")
#for overrides

amplicon_merged_segments_path = paste0(cohort_results_dir,"amplicon_merged_segments.",run_date,".tsv")
amplicon_footprints_path = paste0(cohort_results_dir,"amplicon_footprints.",run_date,".tsv")

ctx_chrom_cn_properties_path = paste0(cohort_results_dir,"ctx_chrom_cn_properties.",run_date,".tsv")
unbalanced_ctx_windows_path = paste0(cohort_results_dir,"unbalanced_ctx_windows.",run_date,".tsv")


## Outputs
#non versioned
complex_sv_bp_cluster_mapping_path = paste0(complex_svs_results_dir,"complex_sv_bp_cluster_mapping.tsv")
#versioned
complex_sv_classification_path = paste0(cohort_results_dir,"complex_sv_classification.",run_date,".tsv")
merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.",run_date,".tsv")
merged_svs_classes_annotated_path = paste0(cohort_results_dir,"merged_svs_classes_annotated.",run_date,".tsv")
complex_svs_annotated_path = paste0(cohort_results_dir,"complex_svs_annotated.",run_date,".tsv")
complex_footprints_annotated_path = paste0(cohort_results_dir,"complex_footprints_annotated.",run_date,".tsv")

map_complex_sv_cn_change_path = paste0(cohort_results_dir,"map_complex_sv_cn_change.",run_date,".tsv")
peaks_complex_sv_cn_change_path = paste0(cohort_results_dir,"recurrent_regions/peaks_complex_sv_cn_change.",run_date,".tsv")
complex_cn_change_regions_path = paste0(cohort_results_dir,"complex_sv_cn_change_regions.",run_date,".tsv")



# Cohort ----

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

automated_baseline_shifts = read.table(automated_baseline_shifts_path,sep="\t",header=T)
manual_baseline_shifts = read.table(manual_baseline_shifts_path,sep="\t",header=T)
baseline_shifts= get_baseline_correction(cohort,automated_baseline_shifts,manual_baseline_shifts)

cancer_type_abbrev = read.table(cancer_type_abbrev_path,sep="\t",header=T)
cohort = cohort %>% left_join(cancer_type_abbrev)


# RESOURCES----

map_template_vars=c('${resources_dir}'=resources_dir)
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

chromosomes_df = get_chromosomes(chr_arms_df)
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames

chr_arm_order = chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")]# %>% flatten_chr()
chrom_order = chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")] %>% unique()# %>% flatten_chr()


# Functions ----
anno_complex_sv_class = function(svs_df,complex_svs_df,manual_labels_path=NULL,
                                 no_complex_class_lst=c("sv_pair","reciprocal_ctx","sv_pair_ctx"),
                                 complex_other_lst=c("complex_other","intra_cluster","ctx_cluster")) {
  #reset class labels, remove cols
  complex_sv_class_labels_cols = c("flag_is_complex","complex_sv_class","complex_sv_class_manual","complex_sv_class_super")
  
  svs_df = svs_df %>% select(names(svs_df)[!names(svs_df) %in% complex_sv_class_labels_cols] )
  complex_svs_df = complex_svs_df %>% select(names(complex_svs_df)[!names(complex_svs_df) %in% complex_sv_class_labels_cols] )
  
  complex_svs_df_cols = names(complex_svs_df)
  complex_svs_df_cols = complex_svs_df_cols[!complex_svs_df_cols %in% c("patient_label","complex_sv_id","cancer_type")]
  
  svs_df = svs_df %>% select(names(svs_df)[!names(svs_df) %in% complex_svs_df_cols] ) %>%  left_join(complex_svs_df)
  svs_df[is.na(svs_df$amplicon_overlap),]$amplicon_overlap=F
  
  if(!is.null(manual_labels_path)) {
    complex_sv_manual_labels = read.table(manual_labels_path,sep = "\t",header = T)
    
    #no complex sv: 
    no_complex_sv = complex_sv_manual_labels[is.na(complex_sv_manual_labels$complex_sv_class),]$complex_sv_id
    
    complex_sv_manual_labels = complex_sv_manual_labels[!is.na(complex_sv_manual_labels$complex_sv_class),]
    #remove if not set
    complex_sv_manual_labels[complex_sv_manual_labels$complex_sv_class=="0",]$complex_sv_class=NA
    complex_sv_manual_labels[complex_sv_manual_labels$remarks=="0",]$remarks=NA
    complex_sv_manual_labels = complex_sv_manual_labels %>% mutate(complex_sv_class_manual=!is.na(complex_sv_class))
    
    svs_df = svs_df %>% left_join(complex_sv_manual_labels,by="complex_sv_id")
    
  } else {
    svs_df$complex_sv_class=NA #make col
  }
  svs_df = svs_df %>% as.data.frame() 
  svs_df = svs_df %>%
    mutate(complex_sv_class = 
             ifelse(!is.na(complex_sv_class),complex_sv_class, #use existing label was based on manual first
                    ifelse(amplicon_overlap,"amplicon", 
                           ifelse(svs_cnt==2, "sv_pair",
                                  ifelse(svs_cnt==3 & ctx_cnt==2 & chrom_cnt==2, "sv_pair_ctx",
                                         ifelse( chrom_cnt==2 & closed_chain & ctx_cnt==4 & dist_chrom_pair_max_mbp<2 & footprint_max_mbp < 20, "reciprocal_ctx",
                                                 ifelse(closed_chain & chrom_cn_balanced_or_edge_all  & !ctx_unbalanced_any,"chromoplexy",
                                                        ifelse(svs_cnt>=10 & svtypes_chisq_pval > 0.05 & seg_cnt_cn_change_sum>=10 & (ctx_chrom_pair_max>=10 | max_sv_cov_complex>=10),
                                                               "chromothripsis",
                                                               "complex_other"))))))))
  
  svs_df$complex_sv_class %>% unique()
  
  #super classes 
  svs_df$complex_sv_class_super=NA
  svs_df = svs_df %>% 
    mutate(complex_sv_class_super = 
             ifelse( (is.na(complex_sv_class) | complex_sv_class %in% no_complex_class_lst),"simple",
                     ifelse(complex_sv_class %in% complex_other_lst,"complex_other",
                            complex_sv_class)))
  
  #override flag is complex based on the class / exclude the simple ones
  svs_df = svs_df %>% mutate(flag_is_complex = ! ( (is.na(complex_sv_class) | complex_sv_class %in% no_complex_class_lst)))
  if(svs_df %>% filter(flag_is_complex & complex_sv_class_super=="simple") %>% nrow() > 0) {
    print("WARNING: ")
    print(svs_df %>% filter(flag_is_complex & complex_sv_class_super=="simple") %>% select(complex_sv_class) %>% unique())
  }
  if(svs_df %>% filter(is.na(complex_sv_class) & complex_sv_class_super!="simple") %>% nrow() > 0) {
    print("WARNING: NAs ")
    
  }
  
  return(svs_df)
}


svtypes_chisq_test = function(target,print_input=F) {
  obs = c(0,0,0,0)
  names(obs) = c("del","inv","inv2","dup")
  obs[1] = sum(target$del_cnt,target$ctx_del_cnt,na.rm = T)
  obs[2] = sum(target$inv_cnt,target$ctx_inv_cnt,na.rm = T)/2
  obs[3] = sum(target$inv_cnt,target$ctx_inv_cnt,na.rm = T)/2
  obs[4] = sum(target$dup_cnt,target$ctx_dup_cnt,na.rm = T)
  chisq = chisq.test(obs, p=rep(1/4,4))
  if(print_input) { print(obs) }
  return(chisq$p.val)
}


# Prepare data ----
## Load SVs ----

unfiltered_svs_df = load_cohort_svs(cohort,c('${merged_svs_dir}'=merged_svs_dir),svs_path_template=svs_union_anno_somatic_multitool_path_template,flag_add_patient_label = T)
unfiltered_svs_df = unfiltered_svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
unfiltered_svs_df = unfiltered_svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
unfiltered_svs_df = annotate_sv_multitool_support(unfiltered_svs_df,sv_tumor_af_threshold = 0.1)
unfiltered_svs_df = unfiltered_svs_df %>% unique()
merged_svs = make_merged_svs(unfiltered_svs_df)
merged_svs = merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  


merged_svs = merged_svs %>% filter(flag_in_filtered_svs) %>% filter(patient_label %in% cohort$patient_label) %>% unique()
#sets globals therefore no output
resolve_orphan_svs(merged_svs,unfiltered_svs_df)

#check consistency merged svs and svs df 
merged_svs %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% nrow() == 0
unfiltered_svs_df %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% filter(patient_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() == 0


merged_svs_gr = GRanges(merged_svs$sv_merged_coordinate)
mcols(merged_svs_gr) = merged_svs
names(merged_svs_gr) = merged_svs_gr$patient_sv_merged


merged_svs_coord = GRanges(merged_svs$sv_merged_coordinate) %>% as.data.frame()
merged_svs_coord$patient_sv_merged = merged_svs$patient_sv_merged
merged_svs = merged_svs %>% left_join(merged_svs_coord[,c("seqnames","start","end","patient_sv_merged")])
merged_svs = merged_svs %>% rowwise() %>% dplyr::mutate(chrom_pair=ifelse(svtype=="CTX",paste0(sort(c(chrom,partner_chrom)),collapse = "-"),NA)) %>% as.data.frame()
merged_svs = merged_svs %>% get_ctx_strands()



## Load CN data ----
segments_df = load_cohort_cn_segments(cohort,map_template_vars=map_template_vars,baseline_shifts=baseline_shifts,flag_add_patient_label=T)
segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id

chrom_centered_segments_df = get_chrom_centered_cn_segments(segments,segments_df,chromosomes,chromosomes_df,cohort) 
chrom_centered_segments = GRanges(chrom_centered_segments_df)
names(chrom_centered_segments) = chrom_centered_segments$patient_cna_id

## Load amplicons ----


amplicon_merged_segments_df = read.table(amplicon_merged_segments_path,sep = "\t",header=T)

amplicon_merged_seg = GRanges(amplicon_merged_segments_df$cna_coordinate)
#mcols(amplicon_merged_seg) = amplicon_merged_segments_df %>% select(-seqnames,-start,-end)
names(amplicon_merged_seg) = amplicon_merged_segments_df$patient_cna_merged_id
amplicon_merged_segments_df = amplicon_merged_segments_df %>% dplyr::rename_with(.cols=c("cr_l2fc_50","footprint_id","cna_width","cna_coordinate"),.fn=function(x){paste0("amplicon_",x)})

amplicon_footprints = read.table(amplicon_footprints_path,sep = "\t",header=T)
amplicon_footprints$footprint_width=amplicon_footprints$footprint_size

amplicon_footprints = amplicon_footprints %>% dplyr::rename_with(.cols=c("footprint_id","footprint_width","footprint_coordinate"),.fn=function(x){paste0("amplicon_",x)})
amplicon_footprints = amplicon_footprints %>% filter(amplicon_footprint_width>amplicons_footprint_size_min)



## CTX balanced/unbalanced ----

ctx_chrom_cn_properties =read.table(ctx_chrom_cn_properties_path,header=T,sep="\t")


# Make complex SV clusters ----

## Load or infer sv clusters ----
# do not override by default

if(length(Sys.glob(complex_sv_bp_cluster_mapping_path))==1 & override_clusters==F) {
  complex_sv_bp_cluster_mapping = read.table(complex_sv_bp_cluster_mapping_path,sep="\t",header = T)
  complex_sv_bp_cluster_mapping = complex_sv_bp_cluster_mapping %>% filter(patient_label %in% cohort$patient_label)
  
} else {
  complex_sv_bp_overlaps = get_complex_sv_bp_overlaps(merged_svs_gr,merged_svs,complex_sv_bp_distance = complex_sv_bp_distance_max) 
  complex_sv_bp_cluster_mapping = get_complex_sv_bp_cluster_mapping(complex_sv_bp_overlaps)
    
  complex_clustered_merged_svs = merged_svs %>% merge(complex_sv_bp_cluster_mapping)
    
  #hashing function used in merge sv but then with merged identifiers instead of sv names
    
  complex_clustered_merged_svs$chrom = factor(complex_clustered_merged_svs$chrom,chrom_order)

  complex_sv_stats = complex_clustered_merged_svs %>% 
      group_by(patient_label,cancer_type,complex_sv_cluster) %>%
      summarize(svs_lst=lst_str(patient_sv_merged),
                svs_cnt=length(unique(patient_sv_merged)),
                partner_svs_lst=lst_str(partner_sv_merged),
                complex_sv_id=paste0("complex_",md5(svs_lst)),
                chrom_lst = toString(sort(unique(chrom))),
                chrom_cnt = length(unique(chrom)),.groups="drop")
    
    complex_sv_stats$complex_sv_id = paste0(complex_sv_stats$patient_label,"_",complex_sv_stats$complex_sv_id)
    
    ## remove those that are just edge ++ partner
    no_complex_sv = complex_sv_stats %>% filter(svs_lst==partner_svs_lst & svs_cnt==2)
    complex_sv_stats = complex_sv_stats %>% filter(!complex_sv_id %in% no_complex_sv$complex_sv_id)
    nrow(complex_sv_stats)
    
    #map back the complex sv id to the cluster mapping and export
    complex_sv_bp_cluster_mapping = complex_sv_bp_cluster_mapping %>% left_join(complex_sv_stats %>% select(-svs_lst,-partner_svs_lst)) %>% filter(!is.na(complex_sv_id))
    
    complex_sv_bp_cluster_mapping = complex_sv_bp_cluster_mapping %>% select(-complex_sv_cluster)
    complex_sv_bp_cluster_mapping$complex_sv_id %>% unique() %>% length()
    
    write.table(complex_sv_bp_cluster_mapping ,complex_sv_bp_cluster_mapping_path,quote=F,sep="\t",col.names = T,row.names = F)
}



## Map back to merged svs ----

merged_svs = merged_svs %>% left_join(complex_sv_bp_cluster_mapping)

## merged grs after anno with complex 
merged_svs_gr = GRanges(merged_svs$sv_merged_coordinate)
mcols(merged_svs_gr) = merged_svs %>% select(-seqnames,-start,-end)
names(merged_svs_gr) = merged_svs_gr$patient_sv_merged

merged_sv_bp_gr = get_svs_as_bp_gr(merged_svs_gr)


## Complex SV local footprints ----

merged_svs$complex_sv_footprint_id = paste0(merged_svs$complex_sv_id,"_",merged_svs$chrom)
merged_svs[is.na(merged_svs$complex_sv_id),]$complex_sv_footprint_id=NA

complex_footprints = merged_svs %>% filter(!is.na(complex_sv_footprint_id)) %>% 
  group_by(patient_label,complex_sv_id,complex_sv_footprint_id) %>% 
  summarize(start=min(start),end=max(end),seqnames=lst_str(seqnames),
            footprint_width=end-start,
            footprint_width_mbp = paste0(round(footprint_width/1e6,2)," Mbp")) %>% 
  get_gr_coordinate(attr_name = "complex_sv_local_footprint_coordinate") %>% ungroup()


complex_footprints_summary = complex_footprints %>% group_by(complex_sv_id) %>%
  summarize(local_footprints = 
              lst_str(paste0(complex_sv_local_footprint_coordinate,"  (",ifelse(footprint_width>1e3,footprint_width_mbp,paste0(footprint_width," bp")),")")),
            footprint_max = max(footprint_width),
            footprint_max_mbp = footprint_max/1e6)


#necessary to use gr? => yes oscillating CN

complex_footprints_gr = GRanges(complex_footprints$complex_sv_local_footprint_coordinate)
names(complex_footprints_gr) = complex_footprints$complex_sv_footprint_id
complex_footprints_gr$complex_sv_footprint_id = complex_footprints$complex_sv_footprint_id



# Infer complex SV features for classification ----

complex_sv_classification = complex_sv_bp_cluster_mapping %>% select(patient_label,cancer_type,complex_sv_id,chrom_cnt,chrom_lst) %>% unique()


## SV types, chromosomes, randomness chisq test ----

# chrom_cnt
# svs_cnt
# ctx_cnt (not adjusted, inter chrom bps)
# svtypes_chisq_pval

complex_sv_stats_svtype_cnts = merged_svs %>% filter(!is.na(complex_sv_id)) %>% 
  group_by(patient_label,complex_sv_id) %>%
  summarize(svs_cnt=length(unique(patient_sv_merged)),
            ctx_cnt=sum(svtype=="CTX"),
            intra_svs_cnt=sum(svtype!="CTX"),
            inv_cnt=sum(svtype=="INV"),
            dup_cnt=sum(svtype=="DUP"),
            del_cnt=sum(svtype=="DEL"))

ctx_types = merged_svs %>% filter(svtype=="CTX") %>% get_ctx_types()
complex_sv_stats_svtype_cnts = complex_sv_stats_svtype_cnts %>% left_join(ctx_types)
complex_sv_stats_svtype_cnts[is.na(complex_sv_stats_svtype_cnts)]=0
complex_sv_stats_svtype_cnts$svtypes_chisq_pval=NA

for(target_id in complex_sv_stats_svtype_cnts$complex_sv_id) {
  test_outcome = complex_sv_stats_svtype_cnts[complex_sv_stats_svtype_cnts$complex_sv_id==target_id,] %>% svtypes_chisq_test() 
  complex_sv_stats_svtype_cnts[complex_sv_stats_svtype_cnts$complex_sv_id==target_id,c("svtypes_chisq_pval")] = test_outcome
}
 
complex_sv_classification = complex_sv_classification %>% left_join(complex_sv_stats_svtype_cnts)

## Closed cycle/chain ----
# closed_chain

complex_sv_bp_overlaps = get_complex_sv_bp_overlaps(merged_svs_gr,merged_svs,complex_sv_bp_distance = complex_sv_bp_distance_max) 
complex_sv_bp_overlaps = annotate_complex_sv_bp_overlaps(complex_sv_bp_overlaps,complex_sv_bp_cluster_mapping,merged_svs)

complex_sv_cycles = annotate_complex_sv_stats_cycles(complex_sv_classification,complex_sv_bp_overlaps,get_cycles_global=T,plot=F)

complex_sv_classification = complex_sv_classification %>% 
  left_join(complex_sv_cycles %>% select(complex_sv_id,contains("cycle"))) %>%
  dplyr::mutate(closed_chain = chrom_cnt>1 & flag_cycle == T & flag_cycle_cnt ==1)


## Distances and densities -----

# dist_chrom_pair_max_mbp

ctx_chrom_pair_gr= merged_sv_bp_gr[merged_sv_bp_gr$svtype=="CTX"]
ctx_chrom_pair_gr$patient_chrom_pair = paste0(ctx_chrom_pair_gr$patient_label,"_",ctx_chrom_pair_gr$chrom_pair)
ctx_chrom_pair_grl = GenomicRanges::split(ctx_chrom_pair_gr,ctx_chrom_pair_gr$patient_chrom_pair)
ctx_chrom_pair_distance_to_nearest_df = get_distance_to_df(ctx_chrom_pair_grl,analyse_sv_bp = T,analysis_type = "nearest")
ctx_chrom_pair_distance_summary = ctx_chrom_pair_distance_to_nearest_df %>% filter(!is.na(complex_sv_id)) %>% group_by(complex_sv_id) %>% get_min_max_median(col="distance_to_nearest",attr_prefix = "dist_chrom_pair_")
ctx_chrom_pair_distance_summary = ctx_chrom_pair_distance_summary %>% select(complex_sv_id,dist_chrom_pair_max) %>% dplyr::mutate(dist_chrom_pair_max_mbp = dist_chrom_pair_max/1e6)

# footprint_max_mbp
# see above complex_footprints_summary

# ctx_chrom_pair_max
ctx_per_chrom_pair = merged_svs %>% filter(svtype=="CTX") %>% group_by(complex_sv_id,chrom_pair) %>% summarize(ctx_cnt=cnt_str(patient_sv_merged))
ctx_per_chrom_pair_summary = ctx_per_chrom_pair %>% group_by(complex_sv_id) %>% get_min_max_median(col="ctx_cnt",attr_prefix="ctx_chrom_pair_")
ctx_per_chrom_pair_summary = ctx_per_chrom_pair_summary %>% select(complex_sv_id,ctx_chrom_pair_max)


# max_cov_complex => max_sv_cov_complex

cov_peaks_per_complex = data.frame()
for(cid in unique(complex_sv_bp_cluster_mapping$complex_sv_id)) {
  complex_gr = merged_svs_gr[!is.na(merged_svs_gr$complex_sv_id) & merged_svs_gr$complex_sv_id==cid]
  cov_peaks = call_peaks(complex_gr,max_peaks_only = F)
  cov_peaks$complex_sv_id=cid
  rownames(cov_peaks)=NULL
  cov_peaks_per_complex = rbind(cov_peaks_per_complex,cov_peaks)
}

max_cov_peaks_per_complex = cov_peaks_per_complex %>% group_by(complex_sv_id) %>% summarize(max_sv_cov_complex=max(cov),total_cov_peaks=n())


complex_sv_classification = complex_sv_classification %>% 
  left_join(ctx_per_chrom_pair_summary) %>%
  left_join(ctx_chrom_pair_distance_summary)  %>% 
  left_join(complex_footprints_summary) %>%
  left_join(max_cov_peaks_per_complex)



## Amplicon overlap ----

# amplicon_overlap

map_amplicons_merged_svs = get_reciprocal_overlap_pairs(resize_gr_distance(amplicon_merged_seg,amplicons_sv_distance_max),merged_svs_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
map_amplicons_merged_svs = map_amplicons_merged_svs %>% dplyr::rename(patient_cna_merged_id=set1,patient_sv_merged=set2) 

#annotate and check same patient
map_amplicons_merged_svs = map_amplicons_merged_svs %>% select(-from,-to) %>%
  left_join(amplicon_merged_segments_df[,c("patient_label","patient_cna_merged_id","amplicon_cna_coordinate","amplicon_cr_l2fc_50","amplicon_footprint_id","amplicon_cna_width")],by="patient_cna_merged_id") %>% 
  left_join(merged_svs[,c("patient_label","patient_sv_merged","complex_sv_id","sv_merged_coordinate")],by="patient_sv_merged")  %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y)

map_amplicons_merged_svs$distance = get_gr_distance(map_amplicons_merged_svs$amplicon_cna_coordinate,map_amplicons_merged_svs$sv_merged_coordinate,flag_as_breakpoints = F) 


#nice to add more info on the amplicon

amplicons_complex_svs_summary = map_amplicons_merged_svs %>% filter(!is.na(complex_sv_id)) %>%
  left_join(amplicon_footprints[,c("amplicon_footprint_id","amplicon_footprint_coordinate","amplicon_footprint_width")],by="amplicon_footprint_id") %>% 
  group_by(complex_sv_id) %>% 
  summarize(#amplicon_cna_merged_lst = lst_str(patient_cna_merged_id),
            amplicon_footprint_lst=lst_str(paste0(amplicon_footprint_id," ",amplicon_footprint_width/1e3," kbp")),
            amplicon_cr_values=lst_str(amplicon_cr_l2fc_50),
            amplicon_cr_mean=mean(amplicon_cr_l2fc_50),
            amplicon_width_kbp=lst_str(paste0(patient_cna_merged_id," ",amplicon_cna_width/1e3," kbp")),
            .groups="drop")


complex_sv_classification = complex_sv_classification %>% 
  dplyr::mutate(amplicon_overlap = complex_sv_id %in% filter(map_amplicons_merged_svs,!is.na(complex_sv_id))$complex_sv_id) %>%
  left_join(amplicons_complex_svs_summary)


## Balanced/unbalanced svs summary ----
# ctx_unbalanced_any
# chrom_cn_balanced_or_edge_all

ctx_chrom_cn_properties = ctx_chrom_cn_properties %>% 
  left_join(complex_sv_bp_cluster_mapping[,c("patient_sv_merged","complex_sv_id")])

ctx_chrom_cn_properties$seqnames = factor(ctx_chrom_cn_properties$seqnames,levels=chrom_order)

ctx_window_properties_complex_summary = ctx_chrom_cn_properties %>% 
  filter(!is.na(complex_sv_id))  %>% 
  group_by(complex_sv_id) %>% #needed to prevent NA group
  summarize(ctx_unbalanced_any = any(chrom_cn_unbalanced),
            ctx_balanced_any = any(chrom_cn_balanced),
            ctx_unbalanced_all = all(chrom_cn_unbalanced),
            ctx_balanced_all = all(chrom_cn_balanced),)

#separate needed is safer
ctx_unbalanced_complex_summary = ctx_chrom_cn_properties %>% 
  filter(!is.na(complex_sv_id) & chrom_cn_unbalanced==T)  %>% 
  group_by(complex_sv_id) %>% #needed to prevent NA group
  summarize(ctx_unbalanced_lst = lst_str(window_feature_cn_value_detailed))

ctx_balanced_complex_summary = ctx_chrom_cn_properties %>% 
  filter(!is.na(complex_sv_id) & chrom_cn_balanced==T)  %>% 
  group_by(complex_sv_id) %>% #needed to prevent NA group
  summarize(ctx_balanced_chrom_lst = lst_str(seqnames)) 

ctx_balanced_or_edge_complex_summary = ctx_chrom_cn_properties %>% 
  filter(!is.na(complex_sv_id) & (chrom_edge | chrom_cn_balanced))  %>% 
  group_by(complex_sv_id) %>% #needed to prevent NA group
  summarize(ctx_balanced_or_edge_chrom_lst = lst_str(seqnames)) 


complex_sv_classification = complex_sv_classification %>% 
  left_join(ctx_window_properties_complex_summary) %>% 
  left_join(ctx_unbalanced_complex_summary) %>%
  left_join(ctx_balanced_complex_summary) %>% 
  left_join(ctx_balanced_or_edge_complex_summary) %>% 
  mutate(chrom_cn_balanced_all=ifelse(!is.na(ctx_balanced_chrom_lst)&ctx_balanced_chrom_lst==chrom_lst,T,F),
         chrom_cn_balanced_or_edge_all=ifelse(!is.na(ctx_balanced_or_edge_chrom_lst)&ctx_balanced_or_edge_chrom_lst==chrom_lst,T,F))


## Oscillating CN ----
# seg_cnt_cn_change_sum


## oscillating pattern:
#get the gains and losses and make mega merge
#then use this to get long format for each gain/loss/0
#because merged know it is not consecutive 

#on combining gain/loss counts across footprints
#choose for summing all gains/loss because that gave the best separation with other sv classes / true positives
#avg doesnt work for multi site ctx like idem ditto max
#sum seems best but maybe over enthousiastic for some complex other that have gains/losses as well.

loss_chrom_centered_segments_merged = chrom_centered_segments[call_cna(chrom_centered_segments)$call=="loss"]
gain_chrom_centered_segments_merged = chrom_centered_segments[call_cna(chrom_centered_segments)$call=="gain"]

loss_chrom_centered_segments_merged = get_cna_merged_gainloss_per_patient(loss_chrom_centered_segments_merged,chrom_centered_segments_df,cnas_df_cols = c("patient_cna_id","cr_l2fc_50","maf_50","call","cna_coordinate","chrom_mean_cr"),group_cols = c("chrom_mean_cr"))
loss_chrom_centered_segments_merged$patient_cna_id=names(loss_chrom_centered_segments_merged)
gain_chrom_centered_segments_merged = get_cna_merged_gainloss_per_patient(gain_chrom_centered_segments_merged,chrom_centered_segments_df,cnas_df_cols = c("patient_cna_id","cr_l2fc_50","maf_50","call","cna_coordinate","chrom_mean_cr"),group_cols = c("chrom_mean_cr"))
gain_chrom_centered_segments_merged$patient_cna_id=names(gain_chrom_centered_segments_merged)

#they shouldnt overlap but just to be sure: 
complex_footprint_cn_chrom_centered_loss = get_chrom_cn_fractions(loss_chrom_centered_segments_merged,loss_chrom_centered_segments_merged %>% as.data.frame(),
                                                                  ranges=complex_footprints_gr,
                                                                  ranges_df=complex_footprints,
                                                                  cohort,return_wide = F,
                                                                  ranges_id_col="complex_sv_footprint_id",ranges_width_col="footprint_width",relative_cr = F)

#complex_footprint_cn_chrom_centered_loss %>% filter(call!="loss")

complex_footprint_cn_chrom_centered_gain = get_chrom_cn_fractions(gain_chrom_centered_segments_merged,gain_chrom_centered_segments_merged %>% as.data.frame(),
                                                                  ranges=complex_footprints_gr,
                                                                  ranges_df=complex_footprints,
                                                                  cohort,return_wide = F,
                                                                  ranges_id_col="complex_sv_footprint_id",ranges_width_col="footprint_width",relative_cr = F)

complex_footprint_cn_chrom_centered = complex_footprints %>% 
  left_join(complex_footprint_cn_chrom_centered_gain %>% pivot_wider(id_cols = c("patient_label","complex_sv_footprint_id"),names_from = "call",values_from = c("seg_cnt","frac_covered"))) %>% 
  left_join(complex_footprint_cn_chrom_centered_loss %>% pivot_wider(id_cols = c("patient_label","complex_sv_footprint_id"),names_from = "call",values_from = c("seg_cnt","frac_covered")) )

complex_footprint_cn_chrom_centered_summary = complex_footprint_cn_chrom_centered %>% group_by(complex_sv_id) %>%
  summarize(
    #seg_cnt_cn_change_avg=mean(c(seg_cnt_loss,seg_cnt_gain),na.rm = T),
    seg_cnt_gain=sum(seg_cnt_gain,na.rm = T),
    seg_cnt_loss=sum(seg_cnt_loss,na.rm = T),
    seg_cnt_cn_change_sum=sum(seg_cnt_loss,seg_cnt_gain,na.rm = T))
    #seg_cnt_cn_change_max=max(seg_cnt_loss,seg_cnt_gain,na.rm = T))

complex_sv_classification = complex_sv_classification %>% 
  left_join(complex_footprint_cn_chrom_centered_summary)
  


# Complex SV classification ----


merged_svs_classes = anno_complex_sv_class(merged_svs,complex_sv_classification,manual_labels_path=manual_labels_path,
                                 no_complex_class_lst=c("sv_pair","reciprocal_ctx","sv_pair_ctx"),
                                 complex_other_lst=c("complex_other","intra_cluster","ctx_cluster"))

complex_sv_class_labels = merged_svs_classes %>% filter(!is.na(complex_sv_id) & !is.na(complex_sv_class)) %>% 
  select(patient_label,cancer_type,flag_is_complex,complex_sv_id,contains("class")) %>% unique() 



complex_sv_classification = complex_sv_classification %>% left_join(complex_sv_class_labels)

complex_footprints = complex_footprints %>% left_join(complex_sv_class_labels)
mcols(complex_footprints_gr)=complex_footprints %>% select(patient_label,complex_sv_id,complex_sv_footprint_id,contains("classes"),flag_is_complex)

# Export ----

write.table(complex_sv_classification,complex_sv_classification_path,sep="\t",col.names = T,row.names = F,quote = T)
write.table(merged_svs_classes,merged_svs_classes_path,sep="\t",col.names = T,row.names = F,quote = T)


# Annotate with genes ----


## Load genes ----

gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)
gene_properties_df = get_gene_properties_df(gtf_path)

gene_properties_df$gene_coord = gene_properties_df$to_coordinate
gene_cols=c("gene_id","ensembl_id","gene_name","gene_type","gene_coord","seqnames")
gene_cols = gene_cols[gene_cols %in% names(gene_properties_df)]


gene_properties=GRanges(gene_properties_df)
names(gene_properties) = gene_properties$gene_id

genes_of_interest_collection = get_genes_of_interest_collection()

gene_properties_df = gene_properties_df %>% 
  mutate(flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
         flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name,
         flag_actionable=gene_name %in% filter(genes_of_interest_collection,db_lst=="actionable_mutations_villani")$gene_name) 

gene_properties_df_flags = names(gene_properties_df)[grepl("flag",names(gene_properties_df))]
for(flag in gene_properties_df_flags) {
  gene_properties_df[is.na(gene_properties_df[,flag]),flag]=F
}


gene_df_cols = c(gene_cols,gene_properties_df_flags)
gene_df_cols[!gene_df_cols %in% names(gene_properties_df)]

genes_1mb = GRanges(gene_properties_df)
names(genes_1mb) =genes_1mb$gene_id
genes_1mb =  resize(genes_1mb, width = width(genes_1mb)+1e6, fix = "center")

genes_names = GRanges(gene_properties_df)
names(genes_names) =genes_names$gene_name



## intermezzo start/end for merged sv > Fusions and especially clinrel fusions ----
##2023-07-31
sv_bp_start_end_gene_overlaps= get_reciprocal_overlap_pairs_start_end(merged_svs_gr,gene_properties,reciprocal_overlap = 0,svtype_matching = FALSE)
sv_bp_start_end_gene_overlaps = sv_bp_start_end_gene_overlaps %>% dplyr::rename(gene_id=set2, patient_sv_merged = set1,svtype=set1_svtype) %>% 
  left_join(gene_properties_df[,gene_df_cols]) %>%  
  left_join(merged_svs_classes %>% dplyr::select(patient_sv_merged,partner_sv_merged,complex_sv_id) %>% unique())

sv_bp_gene_partnered = make_sv_bp_gene_partnered(sv_bp_start_end_gene_overlaps,gene_overlap_cols = gene_df_cols,svs_df_cols = NULL,patient_sv_col="patient_sv_merged",partner_sv_col="partner_sv_merged")

#same with advantage that you also get the CTX partners
sv_bp_gene_wide = sv_bp_gene_partnered %>% filter(flag_gene_of_interest) %>%
  group_by(patient_sv_merged) %>% summarize(start_gene_name = toString(unique(sort(gene_name))),
                                            end_gene_name = toString(unique(sort(partner_gene_name)))) 


sv_bp_gene_wide_pmc = sv_bp_gene_partnered %>% filter(flag_pmc_panel & partner_flag_pmc_panel) %>%
  group_by(patient_sv_merged) %>% summarize(start_gene_name_pmc = toString(unique(sort(gene_name))),
                                            end_gene_name_pmc = toString(unique(sort(partner_gene_name)))) 



## annotate complex svs with genes ----
overlaps_genes_svs_1mb = get_reciprocal_overlap_pairs(genes_1mb,merged_sv_bp_gr,ignore_strand = T,svtype_matching = F,reciprocal_overlap=0) %>%
  dplyr::rename(gene_id=set1, patient_sv_merged = set2 ) %>% 
  left_join(gene_properties_df[,gene_df_cols] %>% dplyr::rename(gene_chrom=seqnames)) %>%  
  left_join(merged_svs_classes %>% dplyr::select(patient_sv_merged,complex_sv_id))


overlaps_genes_svs = get_reciprocal_overlap_pairs(gene_properties,merged_sv_bp_gr,ignore_strand = T,svtype_matching = F,reciprocal_overlap=0) %>%
  dplyr::rename(gene_id=set1, patient_sv_merged = set2 ) %>% 
  left_join(gene_properties_df[,gene_df_cols] %>% dplyr::rename(gene_chrom=seqnames)) %>%  
  left_join(merged_svs_classes %>% dplyr::select(patient_sv_merged,complex_sv_id))


merged_svs_classes_annotated = merged_svs_classes %>%  
  left_join(overlaps_genes_svs %>% filter(flag_gene_of_interest) %>% group_by(patient_sv_merged) %>% summarize(gene_sv_bp_lst=lst_str(gene_name))) %>% 
  left_join(overlaps_genes_svs_1mb %>% filter(flag_gene_of_interest) %>% group_by(patient_sv_merged) %>% summarize(gene_sv_bp_1mb_lst=lst_str(gene_name))) %>%  
  left_join(overlaps_genes_svs %>% filter(flag_pmc_panel) %>% group_by(patient_sv_merged) %>% summarize(gene_sv_bp_pmc_lst=lst_str(gene_name))) %>% 
  left_join(overlaps_genes_svs_1mb %>% filter(flag_pmc_panel) %>% group_by(patient_sv_merged) %>% summarize(gene_sv_bp_1mb_pmc_lst=lst_str(gene_name)))

merged_svs_classes_annotated = merged_svs_classes_annotated %>%
  dplyr::mutate(flag_cancer_gene_sv_bp=!is.na(gene_sv_bp_lst),
                flag_cancer_gene_sv_bp_1mb=!is.na(gene_sv_bp_1mb_lst),
                flag_cancer_gene_sv_bp_pmc = !is.na(gene_sv_bp_pmc_lst),
                flag_cancer_gene_sv_bp_1mb_pmc = !is.na(gene_sv_bp_1mb_pmc_lst),
  )


merged_svs_classes_annotated = merged_svs_classes_annotated %>% left_join(sv_bp_gene_wide) %>% left_join(sv_bp_gene_wide_pmc)
merged_svs_classes_annotated = merged_svs_classes_annotated %>% 
  dplyr::mutate(gene_pair_start_end= ifelse((is.na(start_gene_name) | is.na(end_gene_name) | start_gene_name=="" | end_gene_name=="") |
                                              start_gene_name==end_gene_name,NA,paste0(start_gene_name,"--",end_gene_name)),
                gene_pair_start_end_pmc = ifelse((is.na(start_gene_name_pmc) & is.na(end_gene_name_pmc)) |
                                                   start_gene_name_pmc==end_gene_name_pmc,NA,paste0(start_gene_name_pmc,"--",end_gene_name_pmc)))


write.table(merged_svs_classes_annotated ,merged_svs_classes_annotated_path ,quote=F,sep="\t",col.names = T,row.names = F)

#simplified export view
complex_svs_gene_pairs = merged_svs_classes_annotated %>% filter(!is.na(complex_sv_id)) %>% 
  group_by(complex_sv_id) %>%
  summarize(gene_pair_start_end_lst=lst_str(gene_pair_start_end),
            gene_pair_start_end_pmc_lst=lst_str(gene_pair_start_end_pmc))


complex_svs_annotated = complex_sv_classification %>% filter(flag_is_complex) %>%
  dplyr::select(patient_label,cancer_type,chrom_lst,complex_sv_id,complex_sv_class,complex_sv_class_super) %>%
  left_join(overlaps_genes_svs %>% filter(flag_gene_of_interest) %>% group_by(complex_sv_id) %>% summarize(gene_sv_bp_lst=lst_str(gene_name))) %>% 
  left_join(overlaps_genes_svs_1mb %>% filter(flag_gene_of_interest) %>% group_by(complex_sv_id) %>% summarize(gene_sv_bp_1mb_lst=lst_str(gene_name))) %>%
  left_join(overlaps_genes_svs %>% filter(flag_pmc_panel) %>% group_by(complex_sv_id) %>% summarize(gene_sv_bp_pmc_lst=lst_str(gene_name))) %>% 
  left_join(overlaps_genes_svs_1mb %>% filter(flag_pmc_panel) %>% group_by(complex_sv_id) %>% summarize(gene_sv_bp_1mb_pmc_lst=lst_str(gene_name))) %>%
  left_join(complex_svs_gene_pairs)

write.table(complex_svs_annotated, complex_svs_annotated_path,quote=F,sep="\t",col.names = T,row.names = F)


## annotate footprints with genes ----

overlaps_genes_complex_svs_footprints = get_reciprocal_overlap_pairs(gene_properties,complex_footprints_gr,ignore_strand = T,svtype_matching = F,reciprocal_overlap=0) %>%
  dplyr::rename(gene_id=set1, complex_sv_footprint_id = set2 ) %>%
  left_join(gene_properties_df[,gene_df_cols] %>% dplyr::rename(gene_chrom=seqnames))

complex_svs_local_footprints_cancer_genes = overlaps_genes_complex_svs_footprints %>%
  filter(flag_gene_of_interest) %>% group_by(complex_sv_footprint_id) %>%
  summarize(gene_lst_footprints=lst_str(gene_name))

complex_footprints_annotated = complex_footprints %>%  left_join(complex_svs_local_footprints_cancer_genes)

write.table(complex_footprints_annotated,complex_footprints_annotated_path ,quote=F,sep="\t",col.names = T,row.names = F)


# Get CN changes chrom level  ----
#  version 2023-07-12
## combine ctx and mapped CN changes
# for unbalanced ctxes that are part of complex, get their chrom/arm CN changes with unbalanced_ctx_windows_df
unbalanced_ctx_windows_df =read.table(unbalanced_ctx_windows_path,header=T,sep="\t")
unbalanced_ctx_windows_df = unbalanced_ctx_windows_df %>% left_join(merged_svs_classes %>% select(patient_sv_merged,complex_sv_id,flag_is_complex))

#map to overlapping footprints with map_footprint_gainloss

#2023-06-02 tried with 5mb = works 
mega_merged_seg_gap_size=5e6
#gainloss_chrom_segments_merged_df = get_mega_merged_segments(segments,segments_df,gap_size=5e6,return_as_df = T)

#2023-07-17/18
#using the centered for CN changepoints
#but if I use segments_df then get the actual CN data instead of the shifted so thats better
## use split_by_cntype to get more than just gain/loss
#%>% call_cna_amp_homloss(cr_amplification_threshold = cr_amplification_threshold,cr_hom_loss_threshold = -Inf)
## didnt add value to ue AMP see below

chrom_centered_segments = chrom_centered_segments %>% call_cna()
gainloss_chrom_segments_merged_df = get_mega_merged_segments(chrom_centered_segments,segments_df,gap_size=5e6,return_as_df = T,split_by_cntype = c("gain","loss","loh") )
gainloss_chrom_segments_merged_df = gainloss_chrom_segments_merged_df %>% left_join(cohort[,patient_tumor_id_cols])

#later: rename to get rid of call that can be confusing but this doesnt work!!
#gainloss_chrom_segments_merged_df$patient_cna_id=paste0(gainloss_chrom_segments_merged_df$patient_label,"_mega_merged_",1:length(gainloss_chrom_segments_merged_df))
gainloss_chrom_segments_merged = GRanges(gainloss_chrom_segments_merged_df)
names(gainloss_chrom_segments_merged) = gainloss_chrom_segments_merged$patient_cna_id


#annotate gainloss with chr arms
# dont want to filter just annotate
chr_arms_no_giestain = get_chr_arms(chromosome_bands_df,split_giestain = F) %>% GRanges()
names(chr_arms_no_giestain) = chr_arms_no_giestain$chr_arm

map_gainloss_merged_to_chr_arms = get_reciprocal_overlap_pairs(gainloss_chrom_segments_merged,chr_arms,reciprocal_overlap = 0,svtype_matching = F)
map_gainloss_merged_to_chr_arms = map_gainloss_merged_to_chr_arms %>% dplyr::mutate(patient_cna_id=set1,chr_arm=set2) 
#map_gainloss_merged_to_chr_arms = map_gainloss_merged_to_chr_arms %>% filter(!grepl("acen|gvar|stalk",chr_arm))
gainloss_merged_to_chr_arms = map_gainloss_merged_to_chr_arms %>% group_by(patient_cna_id) %>% summarize(
  overlaps_p_arm=any(grepl("p",chr_arm)),
  overlaps_q_arm=any(grepl("q",chr_arm)),
  chr_arm_lst = lst_str(chr_arm)) #note that  this makes matching/filtering more difficult!

gainloss_chrom_segments_merged_df = gainloss_chrom_segments_merged_df %>% left_join(gainloss_merged_to_chr_arms)
gainloss_chrom_segments_merged_df %>% filter(is.na(overlaps_p_arm) | is.na(overlaps_q_arm)) %>% nrow() == 0

#get complex sv footprint overlaps with mega merged seg
## TODO try actual footprint reduce gr of svs instead of max span on chromosome

map_footprint_gainloss = get_reciprocal_overlap_pairs(complex_footprints_gr,gainloss_chrom_segments_merged,reciprocal_overlap = 0,svtype_matching = F)
map_footprint_gainloss = map_footprint_gainloss %>% dplyr::rename(complex_sv_footprint_id=set1,patient_cna_id=set2) 
map_footprint_gainloss = map_footprint_gainloss %>% dplyr::rename(overlap_complex_cna=overlap_set1_set2,overlap_cna_complex=overlap_set2_set1) 

map_footprint_gainloss = map_footprint_gainloss %>%
  left_join(gainloss_chrom_segments_merged_df %>% 
              select(patient_label,patient_cna_id,call,cna_width,cr_l2fc_50,cna_coordinate,overlaps_p_arm,overlaps_q_arm,chr_arm_lst),by="patient_cna_id") %>%
  left_join(complex_footprints %>% select(patient_label,complex_sv_footprint_id,complex_sv_id),by="complex_sv_footprint_id") %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y) %>%
  left_join(merged_svs_classes %>% select(complex_sv_id,complex_sv_footprint_id,chrom,cancer_type,contains("class"),flag_is_complex) %>% unique() ) %>% unique()

map_footprint_gainloss = map_footprint_gainloss %>% filter(flag_is_complex)

#filtering for good match not accidental overlap
## used min size 5mb but doesnt work for smaller regions eg mycn
#map_footprint_gainloss = map_footprint_gainloss  %>% filter(cna_width>5e6)
# maybe something with min fraction of the footprint overlap relative to the cn change
map_footprint_gainloss = map_footprint_gainloss %>% dplyr::mutate(overlap_cna_bp = overlap_cna_complex*cna_width)
map_footprint_gainloss = map_footprint_gainloss %>% dplyr::rename(cna_mean_cr_l2fc_50=cr_l2fc_50)

#at the moment focus on this, can be others because of relative mega merge seg
map_footprint_gainloss = map_footprint_gainloss %>% filter(call %in% c("gain","loss","loh"))



#build dataframe with simple cols
#cr is of entire window or mega merged seg averaged

map_complex_sv_cn_change = unbalanced_ctx_windows_df %>% 
  filter(call!="0" & flag_is_complex) %>%
  select(patient_label,complex_sv_id,chr_arm,call,window_id,window_coordinate,window_width,window_mean_cr_l2fc_50) %>% 
  dplyr::rename(cna_width=window_width,cna_mean_cr_l2fc_50=window_mean_cr_l2fc_50) %>% unique() %>%
  dplyr::mutate(source="unbalanced_ctx",overlap_cna_bp=NA,overlap_complex_cna=NA,overlap_cna_complex=NA)


#some mega merged seg affect both chr arms 
#add it twice now if overlap both arms... 
chr_arm_level_footprint_gainloss = rbind(map_footprint_gainloss %>% filter(overlaps_p_arm) %>% dplyr::mutate(chr_arm=paste0(chrom,"p")),
                                         map_footprint_gainloss %>% filter(overlaps_q_arm) %>% dplyr::mutate(chr_arm=paste0(chrom,"q")))


chr_arm_level_footprint_gainloss_unfiltered = chr_arm_level_footprint_gainloss
#filter here instead of earlier
chr_arm_level_footprint_gainloss = chr_arm_level_footprint_gainloss_unfiltered  %>% 
  dplyr::mutate(source="footprint_overlap") %>%
  filter(overlap_cna_bp> 5e6 | #substantial overlap  
           #(overlap_complex_cna>0.1 & overlap_cna_complex > 0.1) |
           overlap_cna_complex > 0.5) #most of CN is overlapping 


chr_arm_level_footprint_gainloss = chr_arm_level_footprint_gainloss %>% dplyr::mutate(window_id=patient_cna_id,window_coordinate=cna_coordinate)

#check 
names(chr_arm_level_footprint_gainloss)[names(chr_arm_level_footprint_gainloss) %in% names(map_complex_sv_cn_change)]
names(map_complex_sv_cn_change)[!names(map_complex_sv_cn_change) %in% names(chr_arm_level_footprint_gainloss)]

map_complex_sv_cn_change = rbind(map_complex_sv_cn_change,
                                 chr_arm_level_footprint_gainloss[,names(map_complex_sv_cn_change)] %>% unique())


map_complex_sv_cn_change = map_complex_sv_cn_change %>% left_join(cohort[,c("patient_label","cancer_type")])

write.table(map_complex_sv_cn_change,map_complex_sv_cn_change_path,sep="\t",quote=F,row.names = F,col.names = T)


## get underlying CN data for CN change do to complex SV  ----
complex_svs_annotated = read.table(complex_svs_annotated_path,sep="\t",header=T)

reduced_complex_sv_cn_change  = unlist(GenomicRanges::reduce(split(map_complex_sv_cn_change_gr, map_complex_sv_cn_change_gr$complex_sv_id),ignore.strand=T))
reduced_complex_sv_cn_change$complex_sv_id=names(reduced_complex_sv_cn_change)
reduced_complex_sv_cn_change$range_id = paste0(reduced_complex_sv_cn_change$complex_sv_id,"_reduced_",1:length(reduced_complex_sv_cn_change))
names(reduced_complex_sv_cn_change)=reduced_complex_sv_cn_change$range_id
reduced_complex_sv_cn_change_df = reduced_complex_sv_cn_change %>% as.data.frame()
reduced_complex_sv_cn_change_df = reduced_complex_sv_cn_change_df %>% dplyr::rename(range_width=width)
reduced_complex_sv_cn_change_df = reduced_complex_sv_cn_change_df %>% left_join(complex_svs_annotated %>% select(complex_sv_id,patient_label),by="complex_sv_id")
reduced_complex_sv_cn_change_df = reduced_complex_sv_cn_change_df %>% get_gr_coordinate(attr_name = "range_coordinate")

complex_cn_change_fractions = get_chrom_cn_fractions(segments,segments_df,
                                                     ranges=reduced_complex_sv_cn_change,
                                                     ranges_df=reduced_complex_sv_cn_change_df,
                                                     cohort,return_wide = T,
                                                     ranges_id_col="range_id",ranges_width_col="range_width",relative_cr = F)


complex_cn_change_regions = reduced_complex_sv_cn_change_df %>% select(range_id,complex_sv_id,seqnames,range_coordinate,range_width) %>%
  left_join(complex_svs_annotated %>% select(patient_label,complex_sv_id,cancer_type,complex_sv_class_super)) %>%
  left_join(complex_cn_change_fractions) %>% as.data.frame()

write.table(complex_cn_change_regions,complex_cn_change_regions_path,sep="\t",quote=F,row.names = F,col.names = T)

complex_cn_change_regions  %>%
  filter(cancer_type=="hepatoblastoma" & grepl("chr1",seqnames)) %>% 
  select(complex_sv_id,call,cr_l2fc_50,range_coordinate,range_width) %>% 
  arrange(complex_sv_id,range_coordinate)

#use the cr l2fc values to conclude clonal



## Recurrent regions CN change due to complex SV ----
#Slice and reduce to get overlapping CN change due to complex sv
#2023-07-28
#
#map_complex_sv_cn_change %>% filter(cancer_type=="hepatoblastoma" & chr_arm=="chr1p" & call=="loss")


map_complex_sv_cn_change_gr = GRanges(map_complex_sv_cn_change$window_coordinate)
names(map_complex_sv_cn_change_gr)=map_complex_sv_cn_change$window_id
mcols(map_complex_sv_cn_change_gr)=map_complex_sv_cn_change

#consider splitting per CN change call and per cancer type
split_by_cntype = c("gain","loss","loh")
peaks_complex_sv_cn_change=data.frame()
for(cntype in split_by_cntype) {
  subset=map_complex_sv_cn_change_gr[map_complex_sv_cn_change_gr$call==cntype]
  reduced_complex_sv_cn_change  = unlist(GenomicRanges::reduce(split(subset, subset$patient_label),ignore.strand=T))
  reduced_complex_sv_cn_change$patient_label=names(reduced_complex_sv_cn_change)
  
  peaks = call_peaks(reduced_complex_sv_cn_change,max_peaks_only = F)
  peaks$peak_type="complex_sv_cn_change"
  peaks$call=cntype
  peaks_complex_sv_cn_change = rbind(peaks_complex_sv_cn_change,peaks)
}

peaks_complex_sv_cn_change$peak_id=paste0(peaks_complex_sv_cn_change$peak_id,"_",peaks_complex_sv_cn_change$call)

peaks_complex_sv_cn_change %>% arrange(-cov)
sv_peaks_overlaps = get_sv_peaks_overlaps(svs_gr = map_complex_sv_cn_change_gr,
                                          svs_df = map_complex_sv_cn_change,
                                          peaks_df = peaks_complex_sv_cn_change,
                                          svs_id_col="window_id",peak_id_col="peak_id",
                                          peak_cols = c("peak_id","seqnames","start","end","peak_coordinate","call"),
                                          sv_peak_cols = c("patient_label","window_id","cancer_type","call"))
sv_peaks_overlaps = sv_peaks_overlaps %>% filter(call.x==call.y) #same call type not just bp 
patients_per_peak = get_patients_per_peak(sv_peaks_overlaps)

peaks_complex_sv_cn_change = peaks_complex_sv_cn_change %>% left_join(patients_per_peak)

write.table(peaks_complex_sv_cn_change,peaks_complex_sv_cn_change_path,sep="\t",quote=F,row.names = F,col.names = T)

## look at overlapping regions

peaks_complex_sv_cn_change_gr = GRanges(peaks_complex_sv_cn_change)
peaks_complex_sv_cn_change_gr[is.na(peaks_complex_sv_cn_change_gr$patient_cnt)] %>% length() == 0
peaks_complex_sv_cn_change %>% filter(cov!=patient_cnt) %>% nrow() == 0

peaks_complex_sv_cn_change_gr[grepl("hepatoblastoma",peaks_complex_sv_cn_change_gr$cancer_type_lst) & 
                                peaks_complex_sv_cn_change_gr$call=="loss" &
                                peaks_complex_sv_cn_change_gr$cov>3] %>% reduce()

#with grepl hbl. can also do it per cancer type probably tehn proper region but the approach works!
#chr1:1-34847815 > chr1 1-34880815

peaks_complex_sv_cn_change %>% arrange(-patient_cnt) %>% head()
