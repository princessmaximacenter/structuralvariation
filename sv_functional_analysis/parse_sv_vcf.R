## Merge SVs
## per patient, merge svs from tools to find confident ones
# Tools: manta, delly, gridss

## Last update 2022-02-08
## Note: previously merge_svs.docker.2.R
### now spit up into parse_sv_vcf.somatic.R which collects from 3 tools into ranges file.

#load default.conf 
#optional loading of default.docker.conf 
#load patient conf 
#load unit test conf 

# Config order is important because of overwrites
# general unit test config goes last, which uses patient config for pathnames, 
# also unit test identifier sometimes used in path names 

## Change log 
## Update 2022-02-15
## Refactor, previously merge_svs.docker.2.R
## Changed filtering from tumor AF >0.05; to read support and population db overlap
### now spit up into 
# - parse_sv_vcf.R which collects vcf from 3 tools and prefilters by supporting reads -> sv ranges 
# using sv ranges, first retrieve overlap with population svs 
# - merge_svs.R takes sv ranges and filter out population svs; and merge remaining svs -> sv merged 
# also added filtering log 

## Update ?
## Before merging, first classify SVs in tumor-specific/germline/low af

## update 2021-07-05
# Greatly speed up SVs to ranges with vectorisation
# Intermediate output SVs as ranges
# Single function file: functions.svs.R no longer requires functions from fusion-sq, and no longer uses functions.annotate_sv_genes.R 
# Prefilter VCF to remove low qual, and apply threshold on tumor AF. 

###

## TEST CODE 
if(FALSE){
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.conf');
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.docker.conf');
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.PMCID057AAK.conf');
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf');
  
  #default for sv prefiltering
  read_support_threshold = 6
  
   
}


if(is.null(patient) | is.null(cohort_identifier)) {
  print("EXIT: patient undefined, need config file with tumor/normal ids ")
  quit()
} 

if(!exists("read_support_threshold")) {
  read_support_threshold = 6
}

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(AnnotationDbi)
  library(VariantAnnotation)
  library(StructuralVariantAnnotation)
  library(rtracklayer)
  library(tidyverse)
  library(stringr)
  library(stringi)
})

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.svs.variant_type.R"))

## Fill path templates
## replace data paths here instead of in config to enable overwrites for environments eg local/docker 
map_template_vars=c('${input_vcf_data_dir}'=input_vcf_data_dir,'${merged_svs_dir}'=merged_svs_dir,
                    '${recurrent_output_dir}'=recurrent_output_dir,'${cohort_wdir}'=cohort_wdir,
                    '${cohort_identifier}'=cohort_identifier,'${patient_basename}'=patient$basename)

svs_ranges_path = stri_replace_all_fixed(svs_ranges_path_template,names(map_template_vars), map_template_vars,vectorize=F)

svs_filtering_log_path = stri_replace_all_fixed(svs_filtering_log_path_template,names(map_template_vars), map_template_vars,vectorize=F)

vcf_manta_somatic_path = stri_replace_all_fixed(vcf_manta_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)
vcf_manta_germline_path = stri_replace_all_fixed(vcf_manta_germline_path_template,names(map_template_vars), map_template_vars,vectorize=F)
vcf_delly_path = stri_replace_all_fixed(vcf_delly_path_template,names(map_template_vars), map_template_vars,vectorize=F)
vcf_gridss_path = stri_replace_all_fixed(vcf_gridss_path_template,names(map_template_vars), map_template_vars,vectorize=F)


manta_present=delly_present=gridss_present=T
if(length(Sys.glob(vcf_manta_somatic_path))!=1) {
  print(paste0("Manta somatic missing: ", vcf_manta_somatic_path) )
  manta_present=F
}
if(length(Sys.glob(vcf_manta_germline_path))!=1) {
  print(paste0("Manta diploid missing", vcf_manta_germline_path))
  manta_present=F
}
if(length(Sys.glob(vcf_delly_path))!=1) {
  print(paste0("Delly missing", vcf_delly_path))
  delly_present=F
}
if(length(Sys.glob(vcf_gridss_path))!=1) {
  print(paste0("GRIDSS missing", vcf_gridss_path))
  gridss_present=F
}
#NOTE OK if some are missing


  ## Load vcf
  ## prefilter pass
  
  prefilter_pass = function(vcf_path) {
    prefilterPASS = function(x) { !grepl("MinQUAL|LowQual|LOW_QUAL", x) }
    pre = FilterRules(list(prefilterPASS = prefilterPASS))
    vcf_filtered = filterVcf(vcf_path, "hg38", tempfile(), prefilters=pre)
    return(vcf_filtered)
  }
  
  if(manta_present) vcf_manta_germline_path = prefilter_pass(vcf_manta_germline_path)
  if(manta_present) vcf_manta_somatic_path = prefilter_pass(vcf_manta_somatic_path)
  if(delly_present) vcf_delly_path = prefilter_pass(vcf_delly_path)
  if(gridss_present) vcf_gridss_path = prefilter_pass(vcf_gridss_path)
  
  manta_gr=delly_gr=gridss_gr=GRanges()
  if(manta_present) manta_gr = read_manta_sv_vcf(vcf_manta_germline_path, vcf_manta_somatic_path,patient)
  if(delly_present) delly_gr = read_delly_sv_vcf(vcf_delly_path,patient)
  if(gridss_present) gridss_gr = read_gridss_sv_vcf(vcf_gridss_path,patient)
  
  all_gr = c(manta_gr,delly_gr,gridss_gr)
  all_gr = trim(all_gr)
  
  ## DO NOT filter on unique then you remove the ones that are overlapping between tools because it does not take into account metadata cols
  
  all_gr=all_gr[!duplicated(names(all_gr))]
  
  ## Harmonize
  ## Note: Manta and GRIDSS have  negative svLen for deletions, Delly not. Also Delly can have svlen for CTX? => set to 0 explicitly 
  elementMetadata(all_gr)["svtype"] = get_svtype(all_gr)
  all_gr[all_gr$svtype=="CTX"]$svLen=NA
  all_gr$svLen = abs(all_gr$svLen)
  
  
  #filter sv bp by read support
  metadata = mcols(all_gr) %>% as.data.frame() 
  tumor_var_read_cols= c("tumor_VAR_SR","tumor_VAR_PR","tumor_DV","tumor_RV","tumor_VF")
  tumor_var_read_cols = tumor_var_read_cols[tumor_var_read_cols %in% names(metadata)]
  metadata = metadata %>% mutate(supporting_reads_tumor= rowSums(across(tumor_var_read_cols),na.rm=T)) %>% as.data.frame() 
  
  sv_bp_insufficient_reads = metadata %>% filter(supporting_reads_tumor<(read_support_threshold+1))
  
  mcols(all_gr)=metadata
  all_gr=all_gr[all_gr$supporting_reads_tumor>read_support_threshold]
  
  
  ## OLD: prefilter allelic fraction tumor< 0.05
  #keep manta pass breakpoints without AF too
  #all_gr = all_gr[(is.na(all_gr$tumor_af)|is.na(all_gr$normal_af)&all_gr$FILTER=="PASS"&all_gr$tool=="manta")|
  #                  (!is.na(all_gr$tumor_af)&!is.na(all_gr$normal_af)&all_gr$tumor_af>0.05)]
  
  all_gr = all_gr[all_gr$partner %in% names(all_gr)] #remove unpartnered
  
  ## Make partnered ranges
  
  read_support_colnames = colnames(mcols(all_gr))[grepl("DV|DR|RV|RR|REF|SR|PR|VF",colnames(mcols(all_gr)))]
  read_support_colnames  = read_support_colnames[grepl("tumor|normal",read_support_colnames)]      
  
  sv_metadata_cols =  c("sourceId",  "svtype", "svLen", "partner","FILTER","QUAL", "ALT","REF",
                        "insLen",  "tumor_af", "normal_af",  "somatic", "tool",read_support_colnames)
  
  svs = make_range_svs(all_gr,sv_metadata_cols)
  svs$patient_sv_name = paste0(patient$patient_id,"_",svs$sv_name)
  
  svs_df = as.data.frame(svs)
  
  #Classify variants prior to merging
  svs_df[is.na(svs_df$tumor_af),c("tumor_af")]=0
  svs_df[is.na(svs_df$normal_af),c("normal_af")]=0
  
  svs_df=annotate_sv_af_class(svs_df)
  if(nrow(filter(svs_df,is.na(variant_type)))>0) {
    print("WARNING: some variants could not be classified, this might affect merging")
  }
  
  write.table(svs_df,svs_ranges_path,sep="\t",col.names = T,row.names = F,quote=F)
  
  
  con <- file(svs_filtering_log_path, open="wt")
  writeLines(paste0("gr_delly: ",length(delly_gr)), con)
  writeLines(paste0("gr_manta: ",length(manta_gr)), con)
  writeLines(paste0("gr_gridss: ",length(gridss_gr)), con)
  writeLines(paste0("read_support_threshold: ",read_support_threshold), con)
  writeLines(paste0("sv_bp_insufficient_reads: ",nrow(sv_bp_insufficient_reads)), con)
  writeLines(paste0("svs_after_support_filter: ",nrow(svs_df)), con)
  close(con)
  