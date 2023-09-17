# Driver gene alterations

## Last update 2023-09-03
# make clinrel_mutations_incl_prognostic
# code from solids.analyse_complex_svs.Rmd 
# the output table used to derive complex SV effect and plot oncoprint
# also used as export table
## Notes
# based on manually curated list of driver genes per cancer type 
# note: 2023-07-18 incl prognostic especially for eRMS not clear what exactly drivers are, often multiple mutations
# included poor outcome as well
# definition: advantageous to tumor
# still needs sanity check eg TP53 nearby gain is not disruptive but sv inside is.

## TODO:
# refactor table name to driver_gene_alterations 
# manual exclusions and inclusions better way  with external tables
# export for manuscript including fewer columns

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

#input
map_clinrel_mutations_path = paste0(resources_dir,"map_clinrel_mutations.txt")
run_date=20230701
merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.",run_date,".tsv")

#no version??
patient_gene_overview_detailed_cancer_path = paste0(cohort_results_dir,"gene_level_overview.somatic.detailed.cancer.tsv")

# output
## todo later change to driver_gene_alterations.
driver_gene_alterations_path = paste0(cohort_results_dir,"clinrel_mutations_incl_prognostic.",run_date,".tsv")


# Resources ----
map_clinrel_mutations = read.table(map_clinrel_mutations_path,sep="\t",header=T)
map_clinrel_mutations = map_clinrel_mutations %>% dplyr::mutate(across(where(is.character), trimws))

# Cohort ----
##Load cohort and annotate 

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)

# Read data

## Gene level alteration overview  ----
patient_gene_overview_all_cancer = read.table(patient_gene_overview_detailed_cancer_path,sep="\t",header=T)
patient_gene_overview_all_cancer = patient_gene_overview_all_cancer %>% mutate(alteration_simple=trimws(alteration_simple)) 
patient_gene_overview_all_cancer = patient_gene_overview_all_cancer %>% filter(patient_label %in% cohort$patient_label)
patient_gene_overview_all_cancer = patient_gene_overview_all_cancer %>% dplyr::mutate(patient_gene = paste0(patient_label,"_",gene_name))

#original definition
#patient_gene_overview_detailed = patient_gene_overview_all_cancer %>% filter(flag_affected_by_snv | flag_affected_by_sv | flag_cn_hom_loss | flag_cn_amplification) 

#include 1mb nearby and cn change
driver_gene_alterations = patient_gene_overview_all_cancer %>%
  filter(flag_affected_by_snv | flag_affected_by_sv | flag_cn_hom_loss | flag_cn_amplification |  flag_affected_by_sv_1mb_cn_change) %>%
  merge(map_clinrel_mutations)

#manual_include_driver_gene_alterations = c("M365AAD_MDMD2") no longer needed with 1mb and cn change
manual_exclude_driver_gene_alterations = c("M050AAB_TP53","M050AAB_NF1")

driver_gene_alterations = driver_gene_alterations %>% filter(!patient_gene %in% manual_exclude_driver_gene_alterations)

#manual include
germline_snv_include = patient_gene_overview_all_cancer %>% filter(gene_name=="TP53" & patient_label=="M002AAB")  %>% merge(map_clinrel_mutations) %>% 
  dplyr::mutate(snv_cnt=1,snv="germline",alteration="snv_germline",alteration_simple="snv",flag_affected_by_snv=T,flag_affected=T)

driver_gene_alterations = rbind(driver_gene_alterations,germline_snv_include)

driver_gene_alterations = driver_gene_alterations %>% dplyr::mutate(
  sv_alteration = (flag_affected_by_sv | flag_affected_by_sv_1mb_cn_change), #includes 1mb hom loss and amp too
  complex_sv_alteration = !is.na(complex_sv_id))


## checks ----
#should be 0
driver_gene_alterations %>% filter(complex_sv_1mb_cnt==0 & complex_sv_bp_cnt==0 & complex_sv_alteration) %>% nrow() == 0
driver_gene_alterations %>% filter((complex_sv_1mb_cnt>0 | complex_sv_bp_cnt>0 ) & !complex_sv_alteration) %>% nrow() == 0 
#shouldnt happen anymore now 1mb is included
driver_gene_alterations %>% filter(!sv_alteration & complex_sv_alteration) %>% nrow() == 0
#driver_gene_alterations = driver_gene_alterations %>% dplyr::mutate(complex_sv_alteration=ifelse(!sv_alteration,F,complex_sv_alteration))


# Export ----

if(length(Sys.glob(driver_gene_alterations_path))>0 & override_if_exists==F) {
  print("Output already exists, did not override")
} else {
  write.table(driver_gene_alterations,driver_gene_alterations_path, row.names = F,col.names = T,sep="\t",quote = F)
}

#todo: export for manuscript



