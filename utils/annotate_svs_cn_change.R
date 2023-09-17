# Annotate SVs with CN change around breakpoint of ctxes
## Last update 2023-07-17
## Refactor moved code from elsewhere to here and made two new export tables
# ctx_chrom_cn_properties: 
# - chrom_cn_judgement: final verdict based on chrom_cn_balanced,chrom_edge, chrom_cn_unbalanced. for unbalanced also "before after" notation (ignoring the 0) 
## example values: balanced, edge, gain, loss, gain loss,  inconclusive
# - window_feature_cn_value, easier to parse: based on chrom_or_chr_arm selects the right window_*_feature_value_cn_change
## removed many of the other cols to prevent confusion in this export

##kept feature_value_cn_change and feature_value_cn_balanced in the original ctx properties file for legacy purposes, 
## use window_feature_cn_value instead for most applications
## window_feature_cn_value_detailed to replace feature_value_cn_change but add "inconclusive/balanced/edge"

# unbalanced_ctx_windows_df: for all unbalanced ctx pass, determine chrom/arm and use p/q arm logic to get the terminal segments.
# using get_windows_gr_from_ctx_properties() version always using before or after windows depending on p/q arm

## Try effect of stability threshold 0.7 to 0.66

## Version design 2023-05-05
## Aim: distinguish balanced and unbalanced translocations on the  chromosome (arm) level 
# using a windows based approach, before/after the interchrom breakpoint
## extending to telomeres (chromosome level) or tel/cen (arm level)  => decision
## and double-check with 5mb around the breakpoint

## Description of criteria;
# chromosome CN balanced = windows balanced, and breakpoint balanced or inconclusive,
# chromosome CN unbalanced = windows unbalanced, and breakpoint unbalanced or inconclusive,
# chrom edge = <5mb from teleomere (chrom before after)
# 
# Decisions on state of chromosome/arm/5mb around bp window are the same:
#   CN state of each window depends on highest fraction, then > 0.66 in stable state to be useful information, otherwise inconclusive.
# Difference in call state or >0.2 copy ratio to be different, and less than to be the same (balanced)
# Minimum size for chromosome arms 5mb and chromosome level 10 mb
# 
# TODO: see if it also works if I use the 0.2 always instead of CN state OR copy ratio value difference
# 
# chosing between chromosome level or arm level -> simplified into feature_value_cn_change or feature_value_cn_balanced 
# if balanced/unbalanced detected on chromosome level, choose that one, else choose arm level.


suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
})


wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(wdir,"functions.cna.R"))
source(paste0(wdir,"functions.general.R"))
source(paste0(wdir,"functions.svs.R"))
source(paste0(wdir,"functions.svs.local.R"))

source(paste0(wdir,"solids.conf"))





rundate=20230701 
ctx_window_properties_path = paste0(cohort_results_dir,"ctx_window_properties.",rundate,".tsv")
ctx_chrom_cn_properties_path=paste0(cohort_results_dir,"ctx_chrom_cn_properties.",rundate,".tsv")
unbalanced_ctx_windows_path = paste0(cohort_results_dir,"unbalanced_ctx_windows.",rundate,".tsv")
flag_override = T


## Settings ----
autosomes= c(paste("chr",1:22,sep=""),"chrX")
threshold_window_frac_cn_state=0.66

## Resources ----

map_template_vars=c('${resources_dir}'=resources_dir,'${merged_svs_dir}'=merged_svs_dir,'${utils_output_dir}'=utils_output_dir,
                    '${cohort_wdir}'=cohort_wdir,'${cohort_identifier}'=cohort_identifier)
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband

#2022-09-16: 
# not splitting by giestain to also include acen/gvar/stalk regions for window making 
#2023-01-31: cannot use giestain spit for chromosomes and because of consistency with windows
chr_arms = get_chr_arms(chromosome_bands_df,split_giestain=F)

chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df=as.data.frame(chr_arms)
chr_arms_df = chr_arms_df %>% dplyr::rename(chr_arm_width = width)

## NB: cannot remove the acen/gvar regions from chromosomes because it spans the centromere
#start is always 1 but depends on giestain used or not so this is more flexible
chromosomes_df = chr_arms_df %>% 
  #  filter(!grepl("acen|gvar|stalk|chrY",chr_arm)) %>% 
  group_by(seqnames) %>% summarize(start = min(start),end=max(end),.groups="drop")
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames
chromosomes_df$chrom_width=chromosomes_df$end-chromosomes_df$start+1 #to make equal with chr arm width


## Functions ----

## Read in cohort ----

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]

patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

automated_baseline_shifts = read.table(automated_baseline_shifts_path,sep="\t",header=T)
manual_baseline_shifts = read.table(manual_baseline_shifts_path,sep="\t",header=T)
baseline_shifts= get_baseline_correction(cohort,automated_baseline_shifts,manual_baseline_shifts)
baseline_shifts$flag_is_diploid_no_relative_correction=NA


# Read SVs ----

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


merged_svs_coord = merged_svs_gr %>% as.data.frame()
merged_svs = merged_svs %>% left_join(merged_svs_coord[,c("seqnames","start","end","patient_sv_merged")])
merged_svs = merged_svs %>% rowwise() %>% dplyr::mutate(chrom_pair=ifelse(svtype=="CTX",paste0(sort(c(chrom,partner_chrom)),collapse = "-"),NA)) %>% as.data.frame()
merged_svs = merged_svs %>% get_ctx_strands()


# Read CN segments  ----
segments_df = load_cohort_cn_segments(cohort,map_template_vars=map_template_vars,baseline_shifts=baseline_shifts,flag_add_patient_label=T)
segments_df = segments_df %>% get_gr_coordinate("cna_coordinate")
segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id



# Assess windows chrom and chr arm ----
## just need the windows with the split then calculate frac select highest 
  
#these are the defaults in the function
  windows_df_cols=c("window_id","window_coordinate","window_width")
  window_cn_fractions_cols = c("window_id","call","cr_l2fc_50","seg_covered_cr_stable","cr_stable")
  window_cn_fractions_cols_pivot = window_cn_fractions_cols[!window_cn_fractions_cols %in% c("window_id")]
  #sv_cols=c("patient_label","patient_sv_merged","svtype","svlen_mean","sv_merged_coordinate","tumor_af_mean","partner_chrom","seqnames","chr_arm","partner_chr_arm")
 
  window_annot_cols = c("flag_stable_size_pass","feature_value_before","feature_value_after","flag_cn_change_window","feature_value_cn_change")
  
## make windows first, then assess the windows

##chrom
ctx_windows_chrom_df = get_split_range_by_breakpoint(svs = merged_svs_gr[merged_svs_gr$svtype=="CTX"],
                                                     patients_svs_df = merged_svs %>% filter(svtype=="CTX"),
                                                     ranges=chromosomes,
                                                     ranges_df=chromosomes_df,
                                                     svs_id_col="patient_sv_merged",
                                                     ranges_id_col="seqnames")
ctx_windows_chrom_df = ctx_windows_chrom_df %>% left_join(merged_svs[,c("patient_sv_merged","patient_label")])
  
ctx_chrom_unfiltered =  assess_cn_state_window_before_after(segments,segments_df,
    ctx_windows_chrom_df = ctx_windows_chrom_df,
    svs_id_col="patient_sv_merged",
    patient_id_col="patient_label",
    threshold_window_frac_cn_state = threshold_window_frac_cn_state)

ctx_chrom_unfiltered = ctx_chrom_unfiltered %>% annotate_stable_size_pass_window(unbalanced_ctx_window_width_min = 10e6,min_window_width_both = T)
ctx_chrom_unfiltered = ctx_chrom_unfiltered %>% annotate_cn_change_window(cr_l2fc_threshold_sv_change = 0.2,apply_cr_value_minimum = F)
ctx_chrom_unfiltered = ctx_chrom_unfiltered %>% dplyr::mutate(window_id_before=paste0(window_id_before,"_chrom"),window_id_after=paste0(window_id_after,"_chrom"))
  
  
## chr arm
ctx_windows_chr_arm_df = get_split_range_by_breakpoint(svs = merged_svs_gr[merged_svs_gr$svtype=="CTX"],
                                                     patients_svs_df = merged_svs %>% filter(svtype=="CTX"),
                                                     ranges=chr_arms,
                                                     ranges_df=chr_arms_df,
                                                     svs_id_col="patient_sv_merged",
                                                     ranges_id_col="chr_arm")
ctx_windows_chr_arm_df = ctx_windows_chr_arm_df %>% left_join(merged_svs[,c("patient_sv_merged","patient_label")])

ctx_chr_arm_unfiltered =  assess_cn_state_window_before_after(segments,segments_df,
    ctx_windows_chrom_df = ctx_windows_chr_arm_df,
    svs_id_col="patient_sv_merged",
    patient_id_col="patient_label",
    threshold_window_frac_cn_state = threshold_window_frac_cn_state)
  
ctx_chr_arm_unfiltered = ctx_chr_arm_unfiltered %>% annotate_stable_size_pass_window(unbalanced_ctx_window_width_min = 5e6,min_window_width_both = T)
ctx_chr_arm_unfiltered = ctx_chr_arm_unfiltered %>% annotate_cn_change_window(cr_l2fc_threshold_sv_change = 0.2,apply_cr_value_minimum = F)
ctx_chr_arm_unfiltered = ctx_chr_arm_unfiltered %>% dplyr::mutate(window_id_before=paste0(window_id_before,"_chr_arm"),window_id_after=paste0(window_id_after,"_chr_arm"))
  
 
# Assess region around bp ----

  svs_1mb_windows = make_window_around_sv_bp(svs_gr=merged_svs_gr[merged_svs_gr$svtype=="CTX"],sv_bp_flank_window_size=1e6)


  cn_around_sv_bp =  assess_cn_state_window_before_after(segments,segments_df,
    ctx_windows_chrom_df = svs_1mb_windows,
    svs_id_col="patient_sv_merged",
    patient_id_col="patient_label",
    threshold_window_frac_cn_state = threshold_window_frac_cn_state)
  
  
  #stability but no min size
  cn_around_sv_bp = cn_around_sv_bp %>% annotate_stable_size_pass_window(unbalanced_ctx_window_width_min = 1,min_window_width_both = T)
  cn_around_sv_bp = cn_around_sv_bp %>% annotate_cn_change_window(cr_l2fc_threshold_sv_change = 0.2,apply_cr_value_minimum = F)
  cn_around_sv_bp = cn_around_sv_bp %>% dplyr::mutate(window_id_before=paste0(window_id_before,"_1mb"),window_id_after=paste0(window_id_after,"_1mb"))
  
  cn_around_sv_bp = cn_around_sv_bp %>% 
    dplyr::rename_with(.cols=c(paste0(windows_df_cols,"_before"),paste0(windows_df_cols,"_after")),
                       .fn=function(x){str_replace(x,"window_",paste0("window_1mb_"))}) %>% 
    dplyr::rename_with(.cols=paste0("frac_",c(paste0(window_cn_fractions_cols_pivot,"_before"),paste0(window_cn_fractions_cols_pivot,"_after"))),
                       .fn=function(x){paste0("window_1mb_",x)}) %>% 
    dplyr::rename_with(.cols=all_of(window_annot_cols),
                       .fn=function(x){paste0("window_1mb_",x)}) 
  

  svs_5mb_windows = make_window_around_sv_bp(svs_gr=merged_svs_gr[merged_svs_gr$svtype=="CTX"],sv_bp_flank_window_size=5e6)
  cn_around_sv_bp_5mb =  assess_cn_state_window_before_after(segments,segments_df,
    ctx_windows_chrom_df = svs_5mb_windows,
    svs_id_col="patient_sv_merged",
    patient_id_col="patient_label",
    threshold_window_frac_cn_state = threshold_window_frac_cn_state)
  
  
  #stability but no min size
  cn_around_sv_bp_5mb = cn_around_sv_bp_5mb %>% annotate_stable_size_pass_window(unbalanced_ctx_window_width_min = 1,min_window_width_both = T)
  cn_around_sv_bp_5mb = cn_around_sv_bp_5mb %>% annotate_cn_change_window(cr_l2fc_threshold_sv_change = 0.2,apply_cr_value_minimum = F)
  cn_around_sv_bp_5mb = cn_around_sv_bp_5mb %>% dplyr::mutate(window_id_before=paste0(window_id_before,"_5mb"),window_id_after=paste0(window_id_after,"_5mb"))
  
  cn_around_sv_bp_5mb = cn_around_sv_bp_5mb %>% 
    dplyr::rename_with(.cols=c(paste0(windows_df_cols,"_before"),paste0(windows_df_cols,"_after")),
                       .fn=function(x){str_replace(x,"window_",paste0("window_5mb_"))}) %>% 
    dplyr::rename_with(.cols=paste0("frac_",c(paste0(window_cn_fractions_cols_pivot,"_before"),paste0(window_cn_fractions_cols_pivot,"_after"))),
                       .fn=function(x){paste0("window_5mb_",x)}) %>% 
    dplyr::rename_with(.cols=all_of(window_annot_cols),
                       .fn=function(x){paste0("window_5mb_",x)}) 
  

# CTX window properties - wide df combining chrom and chr arm ----
  
  #characteristics of window -> find-replace window_ by window_chr_arm/chrom
  ## c(paste0(windows_df_cols,"_before"),paste0(windows_df_cols,"_after"))
  
  #values of window => add window_chr_arm/chrom
  ## paste0("frac_",window_cn_fractions_cols_pivot)
  ## paste0("frac_",c(paste0(window_cn_fractions_cols_pivot,"_before"),paste0(window_cn_fractions_cols_pivot,"_after")))
  
  #cols that have been added about feature values => add window_chr_arm/chrom
  
  
  ctx_chr_arm_window_properties = ctx_chr_arm_unfiltered %>% 
    dplyr::rename_with(.cols=c(paste0(windows_df_cols,"_before"),paste0(windows_df_cols,"_after")),
                     .fn=function(x){str_replace(x,"window_",paste0("window_chr_arm_"))}) %>% 
    dplyr::rename_with(.cols=paste0("frac_",c(paste0(window_cn_fractions_cols_pivot,"_before"),paste0(window_cn_fractions_cols_pivot,"_after"))),
                       .fn=function(x){paste0("window_chr_arm_",x)}) %>% 
    dplyr::rename_with(.cols=all_of(window_annot_cols),
                       .fn=function(x){paste0("window_chr_arm_",x)}) 
  
  ctx_chrom_window_properties = ctx_chrom_unfiltered %>% 
    dplyr::rename_with(.cols=c(paste0(windows_df_cols,"_before"),paste0(windows_df_cols,"_after")),
                       .fn=function(x){str_replace(x,"window_",paste0("window_chrom_"))}) %>% 
    dplyr::rename_with(.cols=paste0("frac_",c(paste0(window_cn_fractions_cols_pivot,"_before"),paste0(window_cn_fractions_cols_pivot,"_after"))),
                       .fn=function(x){paste0("window_chrom_",x)}) %>% 
    dplyr::rename_with(.cols=all_of(window_annot_cols),
                       .fn=function(x){paste0("window_chrom_",x)}) 
  
  ctx_window_properties = ctx_chrom_window_properties %>% left_join(ctx_chr_arm_window_properties)
  
  
  ctx_window_properties = ctx_window_properties %>% left_join(cn_around_sv_bp) %>% left_join(cn_around_sv_bp_5mb)
  
  #for interpretation
  ##kept feature_value_cn_change and feature_value_cn_balanced for legacy purposes, use window_feature_cn_value instead for most
  ctx_window_properties =  ctx_window_properties %>% left_join(merged_svs[,c("patient_sv_merged","seqnames","chr_arm")]) %>% 
    dplyr::mutate(feature_value_cn_change = 
                            ifelse( !is.na(window_chrom_flag_cn_change_window) & window_chrom_flag_cn_change_window , 
                                    paste0("chrom ",seqnames,": ",window_chrom_feature_value_cn_change), 
                            ifelse( !is.na(window_chr_arm_flag_cn_change_window) & window_chr_arm_flag_cn_change_window, 
                                    paste0("chr arm ",chr_arm,": ",window_chr_arm_feature_value_cn_change), 
                                    NA)),
                  feature_value_cn_balanced = 
                    ifelse( !is.na(window_chrom_flag_cn_change_window) & !window_chrom_flag_cn_change_window , 
                            paste0("chrom ",seqnames,": ",window_chrom_feature_value_cn_change), 
                            ifelse( !is.na(window_chr_arm_flag_cn_change_window) & !window_chr_arm_flag_cn_change_window, 
                                    paste0("chr arm ",chr_arm,": ",window_chr_arm_feature_value_cn_change), 
                                    NA)))
  
  ctx_window_properties = ctx_window_properties %>%  
    dplyr::mutate(chrom_cn_balanced = !is.na(feature_value_cn_balanced) & (is.na(window_5mb_feature_value_cn_change) | window_5mb_feature_value_cn_change=="balanced"),
                  chrom_cn_unbalanced = !is.na(feature_value_cn_change) & (is.na(window_5mb_feature_value_cn_change) | window_5mb_feature_value_cn_change!="balanced"),
                  chrom_edge = window_chrom_width_before < 5e6 | window_chrom_width_after < 5e6) 
  
  ctx_window_properties[ctx_window_properties$chrom_cn_balanced,c("feature_value_cn_change")]=NA
  ctx_window_properties[ctx_window_properties$chrom_cn_unbalanced,c("feature_value_cn_balanced")]=NA
  
  #2023-07-13 also this to prevent accidental use
  ctx_window_properties[ctx_window_properties$chrom_cn_unbalanced==F,c("feature_value_cn_change")]=NA
  ctx_window_properties[ctx_window_properties$chrom_cn_balanced==F,c("feature_value_cn_balanced")]=NA
  
  #sanity checks:
  #contradicting arms and chrom level should not exist
  ctx_window_properties %>% 
   filter( !is.na(window_chrom_flag_cn_change_window) & window_chrom_flag_cn_change_window ) %>%
   filter( !is.na(window_chr_arm_flag_cn_change_window) & window_chr_arm_flag_cn_change_window==F ) %>% nrow() == 0
  
  ctx_window_properties %>% 
  filter( !is.na(window_chr_arm_flag_cn_change_window) & window_chr_arm_flag_cn_change_window ) %>%
    filter( !is.na(window_chrom_flag_cn_change_window) & window_chrom_flag_cn_change_window==F ) %>% nrow() == 0

  #see overrides above
  ctx_window_properties %>% filter(
    (!is.na(feature_value_cn_change) & !is.na(feature_value_cn_balanced)) | 
    (chrom_cn_balanced & chrom_cn_unbalanced) | 
    (chrom_edge & (chrom_cn_balanced|chrom_cn_unbalanced))
  ) %>% nrow() == 0
  
  #add region selected 
  ctx_window_properties = ctx_window_properties %>%
    dplyr::mutate(chrom_or_chr_arm = ifelse( !is.na(window_chrom_flag_cn_change_window),"chrom",ifelse(!is.na(window_chr_arm_flag_cn_change_window),"chr_arm",NA)))
                

if(flag_override) {
write.table(ctx_window_properties,ctx_window_properties_path,sep="\t",col.names = T,row.names = F,quote = F)
}
  
# Get ctx_chrom_cn_properties, infer unbalanced ctx ----
## rename feature_value_cn_change to something with windows => window_feature_cn_value
## final judgement: chrom_cn_balanced chrom_cn_unbalanced chrom_edge & see what else needed in analysis 

if(any(duplicated(ctx_window_properties$patient_sv_merged))==F & 
   ctx_window_properties %>% filter(!patient_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() == 0 &
   merged_svs %>% filter(svtype=="CTX" & !patient_sv_merged %in% ctx_window_properties$patient_sv_merged) %>% nrow() == 0) {
  print("Pass sanity checks: ctx_window_properties for all merged svs") 
} else {
  warning("inconsistency between ctx_window_properties and merged svs")
}

ctx_chrom_cn_properties = merged_svs %>% filter(svtype=="CTX") %>%
    select(patient_label,patient_sv_merged,partner_sv_merged) %>% 
    left_join(ctx_window_properties, by = c("patient_sv_merged"))
  
  
ctx_chrom_cn_properties = ctx_chrom_cn_properties %>% rowwise() %>%
    dplyr::mutate(window_feature_cn_value = 
                    ifelse( chrom_or_chr_arm=="chrom", window_chrom_feature_value_cn_change,
                            ifelse( chrom_or_chr_arm=="chr_arm", window_chr_arm_feature_value_cn_change, 
                                    NA)),
                  window_cn_call_before  = 
                    ifelse( chrom_or_chr_arm=="chrom", window_chrom_frac_call_before,
                            ifelse( chrom_or_chr_arm=="chr_arm", window_chr_arm_frac_call_before, 
                                    NA)),
                  window_cn_call_after  = 
                    ifelse( chrom_or_chr_arm=="chrom", window_chrom_frac_call_after,
                            ifelse( chrom_or_chr_arm=="chr_arm", window_chr_arm_frac_call_after, 
                                    NA)),
                  chrom_cn_judgement = ifelse(chrom_cn_balanced,"balanced",
                                              ifelse(chrom_edge,"edge",
                                                     ifelse(chrom_cn_unbalanced,
                                                            lst_str(c(trimws(paste(ifelse(window_cn_call_before!="0",window_cn_call_before,""),
                                                                                   ifelse(window_cn_call_after!="0",window_cn_call_after,""),collapse= "-")))),
                                                            "inconclusive"))))

#for display purposes and legacy of feature_value_cn_change
ctx_chrom_cn_properties = ctx_chrom_cn_properties %>%
  dplyr::mutate(window_feature_cn_value_detailed = ifelse(chrom_cn_judgement %in% c("balanced","edge","inconclusive"), chrom_cn_judgement, feature_value_cn_change))
  
ctx_chrom_cn_properties %>% select(chrom_cn_balanced,chrom_cn_unbalanced,chrom_edge,chrom_cn_judgement) %>% unique() %>% 
    arrange(chrom_cn_unbalanced,chrom_cn_balanced,chrom_edge)
  
  
#add partner
ctx_chrom_cn_properties_min_cols = c("chrom_cn_judgement","window_feature_cn_value","chrom_cn_balanced","chrom_cn_unbalanced")

ctx_chrom_cn_properties = ctx_chrom_cn_properties %>% 
    left_join( ctx_chrom_cn_properties[,c("patient_sv_merged",ctx_chrom_cn_properties_min_cols)] %>% 
                 dplyr::rename(partner_sv_merged=patient_sv_merged,
                               partner_chrom_cn_judgement=chrom_cn_judgement,
                               partner_window_feature_cn_value=window_feature_cn_value,
                               partner_chrom_cn_balanced=chrom_cn_balanced,
                               partner_chrom_cn_unbalanced=chrom_cn_unbalanced),
               by=c("partner_sv_merged")) 

ctx_chrom_cn_properties = ctx_chrom_cn_properties %>%  mutate(chr_arm_orientation = ifelse(grepl("p",chr_arm),"p","q"))

export_ctx_chrom_cn_properties_cols= c("patient_label","patient_sv_merged","seqnames","chr_arm","chr_arm_orientation","chrom_or_chr_arm",
                                       ctx_chrom_cn_properties_min_cols,"chrom_edge",
                                       "window_cn_call_before","window_cn_call_after",
                                       "feature_value_cn_change","window_feature_cn_value_detailed", #legacy
                                       "partner_sv_merged",paste0("partner_",ctx_chrom_cn_properties_min_cols))
#check cols removed = OK
names(ctx_chrom_cn_properties)[!names(ctx_chrom_cn_properties) %in% export_ctx_chrom_cn_properties_cols]
export_ctx_chrom_cn_properties_cols[!export_ctx_chrom_cn_properties_cols %in% names(ctx_chrom_cn_properties)]

#subset
export_ctx_chrom_cn_properties_cols = names(ctx_chrom_cn_properties)[names(ctx_chrom_cn_properties) %in% export_ctx_chrom_cn_properties_cols]

if(flag_override) {
write.table(ctx_chrom_cn_properties[,export_ctx_chrom_cn_properties_cols],ctx_chrom_cn_properties_path,sep="\t",col.names = T,row.names = F,quote = F)
}

unbalanced_ctx = ctx_chrom_cn_properties %>% filter(chrom_cn_unbalanced)
  
  
# Unbalanced ctx windows ----
##Get all unb ctx pass, determine chrom/arm 
##For all unb ctx pass, use p/q arm logic to get the terminal segments. => make into GR

sv_property_cols = c("patient_sv_merged")#,"complex_sv_class_super","complex_sv_id","flag_is_complex")

## maybe add frac_cr_l2fc_50 from the unbalanced ctx
unbalanced_ctx_gr_chrom = get_windows_gr_from_ctx_properties(unbalanced_ctx %>% filter(chrom_or_chr_arm=="chrom"),region = "chrom",property_cols = sv_property_cols)
unbalanced_ctx_gr_chr_arm = get_windows_gr_from_ctx_properties(unbalanced_ctx %>% filter(chrom_or_chr_arm=="chr_arm"),region = "chr_arm",property_cols = sv_property_cols)
  
unbalanced_ctx_gr = c(unbalanced_ctx_gr_chrom,unbalanced_ctx_gr_chr_arm)
unbalanced_ctx_windows_df = unbalanced_ctx_gr %>% as.data.frame()

unbalanced_ctx_windows_df$chr_arm = paste0(unbalanced_ctx_windows_df$seqnames,unbalanced_ctx_windows_df$chr_arm_orientation)
unbalanced_ctx_windows_df = unbalanced_ctx_windows_df %>% get_gr_coordinate(attr_name = "window_coordinate")
unbalanced_ctx_windows_df = unbalanced_ctx_windows_df %>% dplyr::rename(window_width=width)

#add avg cn over the windows
window_chrom_cn_fractions = get_chrom_cn_fractions(segments,segments_df, 
                                                   ranges=unbalanced_ctx_gr,
                                                   ranges_df=unbalanced_ctx_windows_df,
                                                   cohort,return_wide = T,
                                                   ranges_id_col="window_id",ranges_width_col="window_width",relative_cr = F)


unbalanced_ctx_windows_df = unbalanced_ctx_windows_df %>% left_join(
  window_chrom_cn_fractions %>% dplyr::rename(window_mean_cr_l2fc_50=cr_l2fc_50) %>% select(window_id,window_mean_cr_l2fc_50) )


if(flag_override) {
  write.table(unbalanced_ctx_windows_df,unbalanced_ctx_windows_path,sep="\t",col.names = T,row.names = F,quote = F)
}

# check missed examples ----

ctx_chrom_cn_properties  %>% filter(patient_label %in% c("M260AAB", "M454AAD", "M909AAA")) %>% filter(seqnames=="chr17") %>%
  select(-names(ctx_window_properties))

unbalanced_ctx_windows_df  %>% filter(patient_label %in% c("M260AAB", "M454AAD", "M909AAA")) %>% filter(seqnames=="chr17")
#M260AAB_merged_a9568cbe6018ed16b33f9b75cbf0e164 before window 68% stable instead of 0.7
#same for M909AAA_merged_dd0428dcc8750e51bb355951ca553170
#M454AAD_merged_c2ae948ce08e7d3aad7e93036ff23dc0 stable after 0.47 but also loss?? 
#=> yes that is correct, I do think the 'after' neutral should pass tho on arm level. but becuase bp is is in q arm it doesnt check for p loss only 
ctx_chrom_cn_properties  %>% filter(patient_sv_merged=="M454AAD_merged_c2ae948ce08e7d3aad7e93036ff23dc0") %>% as.data.frame()


ctx_chrom_cn_properties %>% 
  filter( (window_chrom_frac_cr_stable_before & window_chrom_frac_cr_stable_after==F & window_chrom_frac_seg_covered_cr_stable_after >0.65) | 
            (window_chrom_frac_cr_stable_after & window_chrom_frac_cr_stable_before==F & window_chrom_frac_seg_covered_cr_stable_before>0.65) | 
           ( window_chr_arm_frac_cr_stable_before & window_chr_arm_frac_cr_stable_after==F & window_chr_arm_frac_seg_covered_cr_stable_after >0.65) | 
  (window_chr_arm_frac_cr_stable_after & window_chr_arm_frac_cr_stable_before==F & window_chr_arm_frac_seg_covered_cr_stable_before>0.65)) %>%
  filter(window_5mb_flag_cn_change_window) %>% select(chrom_cn_judgement)

# Test cases to check values ----
  if(FALSE){
  new_ctx_window_properties =read.table(paste0("~/results/solids_v4/ctx_window_properties.20230505.tsv"),header=T,sep="\t")
  merged_svs_classes = read.table(paste0(cohort_results_dir,"merged_svs_classes.20230516.tsv"),sep="\t",header=T)
  
  #differences:
  ctx_window_properties %>% anti_join(new_ctx_window_properties[,c("patient_sv_merged","feature_value_cn_change","feature_value_cn_balanced")])
  ctx_window_properties %>% anti_join(new_ctx_window_properties[,c("patient_sv_merged","chrom_cn_balanced","chrom_cn_unbalanced")])
  
  diff = ctx_window_properties %>% left_join(new_ctx_window_properties[,c("patient_sv_merged","feature_value_cn_change","feature_value_cn_balanced",
                                                                   "chrom_cn_balanced","chrom_cn_unbalanced")],by="patient_sv_merged")
  
  diff %>% filter(chrom_cn_balanced.x!=chrom_cn_balanced.y |
                    chrom_cn_unbalanced.x!=chrom_cn_unbalanced.y) %>% 
    #filter(chrom_cn_balanced.y==F & chrom_cn_unbalanced.y==F) %>% #these were inconclusive
    filter(chrom_cn_balanced.y==T | chrom_cn_unbalanced.y==T) %>% 
    select(patient_sv_merged,contains("feature_value_cn_"),contains("balanced")) %>% 
    left_join(merged_svs_classes %>% select(patient_sv_merged,contains("complex"))) %>% View()
  
  
  #some became inconlusive others became conclusive => effect on complex svs? = > No
  }


if(FALSE){
  prev_ctx_chrom_cn_properties =read.table(paste0("~/results/solids_v4/ctx_chrom_cn_properties.20230701.tsv"),header=T,sep="\t")
  merged_svs_classes = read.table("~/results/solids_v4/merged_svs_classes.20230701.tsv",sep="\t",header=T)
  
  #differences:
  ctx_chrom_cn_properties = ctx_chrom_cn_properties %>% select(export_ctx_chrom_cn_properties_cols) %>% anti_join(prev_ctx_chrom_cn_properties[,c("patient_sv_merged","chrom_cn_judgement")])
  ctx_chrom_cn_properties %>% select(export_ctx_chrom_cn_properties_cols) %>% anti_join(prev_ctx_chrom_cn_properties[,c("patient_sv_merged","chrom_cn_balanced","chrom_cn_unbalanced")])
  
  diff = ctx_chrom_cn_properties %>% select(export_ctx_chrom_cn_properties_cols) %>% 
    left_join(prev_ctx_chrom_cn_properties[,c("patient_sv_merged","chrom_cn_judgement","chrom_cn_balanced","chrom_cn_unbalanced")],by="patient_sv_merged")
  
  diff %>% filter(chrom_cn_balanced.x!=chrom_cn_balanced.y |
                    chrom_cn_unbalanced.x!=chrom_cn_unbalanced.y) %>% 
    filter(chrom_cn_balanced.y==F & chrom_cn_unbalanced.y==F) %>% #these were inconclusive
    #filter(chrom_cn_balanced.y==T | chrom_cn_unbalanced.y==T) %>% 
    select(patient_sv_merged,seqnames, window_feature_cn_value,contains("balanced")) %>% 
    left_join(merged_svs_classes %>% select(patient_sv_merged,contains("complex"))) %>% View()
}
  
  #ctx_window_properties %>% filter(patient_sv_merged=="M974AAC_merged_6ac12cc58362f46463f41f672138a570")

  #test case: cn_around_sv_bp_5mb %>% filter(patient_sv_merged == "M228AAA_merged_aface363c6bddf7b887f3a8ca9487559")
  #should be neither balanced nor balanced because windows tell another story then around breakpoint
  #ctx_window_properties %>% filter(patient_sv_merged=="M228AAA_merged_aface363c6bddf7b887f3a8ca9487559")

  #test case: not balanced and stable, nearby clear arm level change
  #is 1mb too small? M075AAD_complex_f0cdf60d48fce595ffc1baaf5bd6cacc Id say still unbalanced, just another inv nearby
  #cn_around_sv_bp %>% filter(patient_sv_merged == "M075AAD_merged_023fccc1ff768ac7569a4bcb395a2eea")
  #cn_around_sv_bp_5mb %>% filter(patient_sv_merged == "M075AAD_merged_023fccc1ff768ac7569a4bcb395a2eea")
  ## should NOT be balanced: M228AAA_complex_70318304349a7ae1451e4884a1e3b975 chr 4 M228AAA_merged_aface363c6bddf7b887f3a8ca9487559

  #edge  => <5mb of chrom before or after means can never pass chrom or arm level
  #M226AAD_merged_2c638659dac11ecd48aa6d8a44850968
  #M226AAD_merged_67e5194d825133602e095005721fd23f
  
