# Complex SVs effect

## Last update 2023-09-03

## TODO:
# manual exclusions and inclusions better way with external tables
# export for manuscript including fewer columns & additional properties
# refactor clinrel CN change to elsewhere? like driver gene alterations /clinrel mutations 

wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(script_dir,"functions.general.R"))

suppressPackageStartupMessages({
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
})

source(paste0(wdir,"solids.conf"))


override_if_exists = F
minimum_cna_width = 5e6


#input
map_clinrel_chromalt_path = paste0(resources_dir,"map_clinrel_chromalt.txt")

run_date=20230701
merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.",run_date,".tsv")
merged_svs_classes_annotated_path = paste0(cohort_results_dir,"merged_svs_classes_annotated.",run_date,".tsv")
complex_svs_annotated_path = paste0(cohort_results_dir,"complex_svs_annotated.",run_date,".tsv")
complex_sv_classification_path = paste0(cohort_results_dir,"complex_sv_classification.",run_date,".tsv")

map_complex_sv_cn_change_path = paste0(cohort_results_dir,"map_complex_sv_cn_change.",run_date,".tsv")

## todo later change to driver_gene_alterations.
driver_gene_alterations_path = paste0(cohort_results_dir,"clinrel_mutations_incl_prognostic.",run_date,".tsv")

# output
curated_complex_sv_cn_change_path = paste0(cohort_results_dir,"map_complex_sv_cn_change.",run_date,".curated.tsv")

complex_svs_effect_path = paste0(cohort_results_dir,"complex_svs_annotated.",run_date,".effect.tsv")

# Resources ----

## Chromosomes ----
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)

chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)
chr_arms_df$chr_arm_width=chr_arms_df$width

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband


chromosomes_df = get_chromosomes(chromosome_bands_df)
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames

chr_arm_order = chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")]# %>% flatten_chr()
chrom_order = chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")] %>% unique()# %>% flatten_chr()


## target list ----
map_clinrel_chromalt = read.table(map_clinrel_chromalt_path,sep="\t",header=T)
map_clinrel_chromalt = map_clinrel_chromalt %>% dplyr::mutate(across(where(is.character), trimws))
map_clinrel_chromalt = map_clinrel_chromalt %>% mutate(chr_arm_cn_change=paste(chr_arm,cn_change,sep = "_"))

# Cohort ----
##Load cohort and annotate 

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)

# Read data ----

## complex svs ----

merged_svs = read.table(merged_svs_classes_path,sep="\t",header=T)
complex_svs_annotated = read.table(complex_svs_annotated_path,sep="\t",header=T)
complex_sv_classification = read.table(complex_sv_classification_path,sep="\t",header=T)

## driver gene alterations ----
driver_gene_alterations = read.table(driver_gene_alterations_path,sep="\t",header=T)

driver_gene_alteration_complex_summary= driver_gene_alterations %>% filter(complex_sv_alteration) %>% 
  group_by(patient_label,complex_sv_id) %>% summarize(driver_gene_lst=lst_str(gene_name),sv_driver_expected=any(sv_expected),complex_sv_driver_expected=any(complex_sv_expected))

if("driver_gene_lst" %in% names(complex_svs_annotated)) {
  complex_svs_annotated = complex_svs_annotated %>% select(-driver_gene_lst,-sv_driver_expected,-complex_sv_driver_expected)
}
complex_svs_annotated = complex_svs_annotated %>% left_join(driver_gene_alteration_complex_summary)

complex_svs_annotated[is.na(complex_svs_annotated$sv_driver_expected),c("sv_driver_expected")]=F
complex_svs_annotated[is.na(complex_svs_annotated$complex_sv_driver_expected),c("complex_sv_driver_expected")]=F


## clinrel CN change -----
map_complex_sv_cn_change = read.table(map_complex_sv_cn_change_path,sep="\t",header = T)
map_complex_sv_cn_change = map_complex_sv_cn_change %>% left_join(complex_svs_annotated %>% select(complex_sv_id,complex_sv_class) %>% unique())

map_complex_sv_cn_change = map_complex_sv_cn_change %>%
  mutate(chr_arm_cn_change=paste(chr_arm,call,sep = "_"),
         patient_chr_alt=paste(patient_label,chr_arm,call,sep = "_"))

map_complex_sv_cn_change$chr_arm = factor(map_complex_sv_cn_change$chr_arm,levels=chr_arm_order)

#2023-08-17 remove manual have seen not LOH due to complex but large footprint overlap and is amplicon type
manual_no_complex_cn_change_lst = c("M002AAB_complex_249866b9973feb77be459324d1531e38")
map_complex_sv_cn_change = map_complex_sv_cn_change %>% filter(!complex_sv_id %in% manual_no_complex_cn_change_lst)

#remove too small to be relevant
map_complex_sv_cn_change = map_complex_sv_cn_change %>% filter(cna_width>minimum_cna_width)

clinrel_complex_sv_cn_change = map_complex_sv_cn_change %>% merge(map_clinrel_chromalt) %>% filter(call==cn_change)

map_complex_sv_cn_change = map_complex_sv_cn_change %>% mutate(clinrel=patient_chr_alt %in% clinrel_complex_sv_cn_change$patient_chr_alt)

#beware of CN change due to amplicon, have seen a few cases where it is just underlying aneuploidy, but can still be from the rearrangement so cant filter out all 
map_complex_sv_cn_change %>% filter(complex_sv_class=="amplicon") %>% select(complex_sv_id,chr_arm,call,source,clinrel) %>% unique() %>% arrange(complex_sv_id)

if(FALSE) {
  #to explore:
  map_complex_sv_cn_change = map_complex_sv_cn_change %>% 
  dplyr::mutate(cn_details_display = paste0(chr_arm," ",call," (",
ifelse(source=="footprint_overlap",paste0(round(overlap_cna_bp/1e6,2)," Mbp of "),
       "unb ctx "),round(cna_width/1e6,2)," Mbp @ ",round(cna_mean_cr_l2fc_50,2)," crl2fc)"))

#lot of duplication, do not include
complex_sv_footprint_cn_change_summary = map_complex_sv_cn_change %>% 
  filter(source=="footprint_overlap") %>% group_by(patient_label,complex_sv_id) %>% 
  summarize(footprint_overlap_cn_change = lst_str(paste0(chr_arm," ",call," (",cn_details_display,")")))

#use the reduced version maybe? but then ranges and not the separate segs
complex_cn_change_regions_path = paste0(cohort_results_dir,"complex_sv_cn_change_regions.",run_date,".tsv")
complex_cn_change_regions = read.table(complex_cn_change_regions_path,sep="\t",header=T)

}

##todo rename to clinrel_cn_change
clinrel_complex_sv_cn_change_summary = clinrel_complex_sv_cn_change %>% 
  group_by(patient_label,complex_sv_id) %>% 
  summarize(driver_cn_change = lst_str(paste0(chr_arm," ",cn_change," (",source,")")))

if("driver_cn_change" %in% names(complex_svs_annotated)) {
  complex_svs_annotated = complex_svs_annotated %>% select(-driver_cn_change) 
}
complex_svs_annotated = complex_svs_annotated %>% left_join(clinrel_complex_sv_cn_change_summary)

## Infer complex effect ----
complex_svs_annotated = complex_svs_annotated %>% dplyr::mutate(
  effect = ifelse(!is.na(driver_gene_lst)&!is.na(driver_cn_change),"driver_and_chrom_alt",
            ifelse(!is.na(driver_gene_lst),"driver",
            ifelse(!is.na(driver_cn_change),"chrom_alt","unknown"))))


## add for manuscript supplemental ----
#properties from classification, tumor allele fraction, 

tumor_af_complex = merged_svs %>% filter(flag_is_complex) %>% group_by(patient_label,complex_sv_id) %>%
  summarize(sv_tumor_af_mean=mean(tumor_af_mean,na.rm=T),sv_tumor_af_max=max(tumor_af_mean,na.rm = T))

complex_svs_annotated = complex_svs_annotated %>% left_join(tumor_af_complex)

if(FALSE) {
  #minimum
  #but want all properties for classification 
selected_cols = c("amplicon_overlap","svs_cnt","ctx_cnt","chrom_cnt",
                  "closed_chain","dist_chrom_pair_max_mbp","footprint_max_mbp",
                  "chrom_cn_balanced_or_edge_all","ctx_unbalanced_any",
                  "svtypes_chisq_pval","seg_cnt_cn_change_sum","ctx_chrom_pair_max","max_sv_cov_complex")
}

                                       
# Export ----
if(length(Sys.glob(complex_svs_effect_path))>0 & override_if_exists==F) {
  print("Output already exists, did not override")
} else {
  write.table(complex_svs_annotated,complex_svs_effect_path,row.names = F,col.names = T,sep="\t",quote = F)
    
  #write for oncoprint, because manually curated
  #todo refactor later?
  write.table(map_complex_sv_cn_change,curated_complex_sv_cn_change_path, row.names = F,col.names = T,sep="\t",quote = F)
  
  
}

#todo: export for manuscript



