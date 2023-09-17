## Get suggestion for baseline shift correction based on stable chr arms 
## goal is to have 2n as copy ratio 0
## suggestion is based on lowest cr stable maf stable chromosome arm
## also return mean and median of chr called as neutral 

## Last update: 2022-10-12


if(FALSE) {
  #local
  wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
  source(paste0(wdir,"default.conf"))
  source(paste0(wdir,"solids.conf"))
  
  run_biomaterial_id="PMLBM000EFX"
  run_patient_id="PMCID576AAQ"
}

if(FALSE){
  #config example 
  
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.conf")
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.docker.conf")
   
  run_biomaterial_id="PMABM000GFX" #required
  run_patient_id="PMCID270AAA"   #add patient id but do not rely on it being there?.. not sure if sensible
  segments_path = "" #sys glob based on biomaterial id if not provided
}


#sys glob with biomaterial id if no path provided 
if(!exists("run_biomaterial_id")) {
  print("No biomaterial id provided, exiting")
  quit()
}


source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.cna.R"))


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

##

map_template_vars=c('${resources_dir}'=resources_dir,'${utils_output_dir}'=utils_output_dir,'${biomaterial_id}'=run_biomaterial_id)


#output paths ----
chr_arm_cna_path = stri_replace_all_fixed(chr_arm_cna_path_template,names(map_template_vars), map_template_vars,vectorize=F)

cna_baseline_shift_path_template = "${utils_output_dir}cna_baseline_shift.${biomaterial_id}.tsv"
cna_baseline_shift_path = stri_replace_all_fixed(cna_baseline_shift_path_template,names(map_template_vars), map_template_vars,vectorize=F)

if(length(Sys.glob(chr_arm_cna_path))==0) {
  print("Input not found:")
  print(chr_arm_cna_path)
  quit()
}
  
#Read CN chr arm data ---- 

chr_arm_cna = read.table(chr_arm_cna_path,sep="\t",header=T)
if(FALSE) {
  ## 2023 06 29 todo finish this
cohort = c()
cohort$patient_id=run_patient_id
cohort$tumor_id=run_biomaterial_id
cohort$patient_label=cohort$patient_id
cohort=as.data.frame(cohort)
segments_df = load_cohort_cn_segments(cohort,map_template_vars=map_template_vars,baseline_shifts=NULL,flag_add_patient_label=T)
segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id

# Chromosomes
map_template_vars=c('${resources_dir}'=resources_dir)
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)
chr_arms_df = get_chr_arms(chromosome_bands_df)
chromosomes_df = get_chromosomes(chr_arms_df)
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames

chr_arm_cna = get_chrom_cn_fractions(segments,segments_df,
                                            chromosomes,chromosomes_df,cohort,ranges_id_col="seqnames",ranges_width_col="chrom_width",
                                            relative_cr = F,return_wide = T)
}
#3 copies (1.5x = 0.58 l2fc, 4n=2x=1 l2fc)
threshold_maf_hom=0.47
min_cr_range_threshold=0.33
#suggestion for correction of shift
#- select chr arms that are cr and maf stable and MAF normal.
stable_chrom =  chr_arm_cna %>% filter(!grepl("_",chr_arm)&!grepl("chrX|chrY",chr_arm)) %>% filter(cr_stable & maf_stable) %>% filter(maf_50>threshold_maf_hom) 

#- select lowest, make cr window +- 10% of that value
lowest_cr_stable_chrom = min(stable_chrom$cr_l2fc_50)
min_cr_range = ifelse(lowest_cr_stable_chrom<0,-min_cr_range_threshold,min_cr_range_threshold) #switch the sign
min_cr = (1-min_cr_range)*lowest_cr_stable_chrom
max_cr = (1+min_cr_range)*lowest_cr_stable_chrom

stable_chrom_lowest = stable_chrom %>% filter(  cr_l2fc_50 < max_cr &  cr_l2fc_50 > min_cr)
#- mean copy ratio of of chr arms (maf and cr stable) in that range
#- suggested correction is flip the sign of observed shift > if the mean is below zero than add the value to correct

suggestion_chrom_summary = stable_chrom_lowest %>% group_by(patient_id) %>% summarize(
  baseline_shift_observed=mean(cr_l2fc_50,na.rm=T),
  shift_chrom_lst=toString(paste0(chr_arm," (",round(cr_l2fc_50,2)," ",round(maf_50,2),")" ))
)  

# add the mean/median copy ratio of neutral chrom
# if centered around 0 then indication that shift is not needed
neutral_stable_chrom_summary = stable_chrom %>% filter(call=="0") %>% group_by(patient_id) %>% summarize(
  mean_neutral_cr=mean(cr_l2fc_50,na.rm=T),
  median_neutral_cr=median(cr_l2fc_50,na.rm=T),
  neutral_chrom_lst=toString(paste0(chr_arm," (",round(cr_l2fc_50,2)," ",round(maf_50,2),")" ))
  )  


cna_baseline_shift = suggestion_chrom_summary %>% left_join(neutral_stable_chrom_summary)
cna_baseline_shift$biomaterial_id=run_biomaterial_id

write.table(cna_baseline_shift,cna_baseline_shift_path,sep="\t",quote=F,row.names = F, col.names = T)
  