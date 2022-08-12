## Pairwise overlap 

# Config is loaded upfront 

if(FALSE){
  #using config or paths
  
  source("~/structuralvariation/sv_functional_analysis/default.conf")
  source("~/structuralvariation/sv_functional_analysis/default.docker.conf")
  source("~/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf")
  source("~/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.patient.conf")
  
  #alternatively specify 
  #svs_path=""
  #patient_id=""
  
}

suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
})
source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))

## Fill path templates
map_template_vars=c('${resources_dir}'=resources_dir,'${merged_svs_dir}'=merged_svs_dir,'${utils_output_dir}'=utils_output_dir,
                    '${cohort_wdir}'=cohort_wdir,'${cohort_identifier}'=cohort_identifier,'${patient_basename}'=patient$basename)

svs_union_path = stri_replace_all_fixed(svs_union_path_template,names(map_template_vars), map_template_vars,vectorize=F)
sv_exon_overlaps_path = stri_replace_all_fixed(sv_exon_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)



if(!exists("svs_path")){
  svs_path=svs_union_path
}
if(!exists("patient_id")) {
  if(exists("patient")){
    patient_id=patient$patient_id
  } else {
    print("No patient id, or patient$ is specified")
    quit()
  }
}


## loading can be standardized across utils scripts
svs_df = read_svs_df(svs_path,patient_id=patient_id)
svs_gr = get_svs_gr(svs_df)

#svs_df_overlap_cols defined in default.conf but can be overwritten
#include for interpretation
#removed patient_sv_merged of reporting because it can differ with version updates 
#svs_df_overlap_cols = c("patient_sv_name","from_coordinate","svLen",
#                "sv_name","bp_names_origin","patient_id","partner","tool","tumor_af","normal_af"))

#as fail safe, only report columns that exist
svs_df_overlap_cols = svs_df_overlap_cols[svs_df_overlap_cols %in% names(svs_df)]

### RESOURCES 
gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)

exon_properties_df = get_exon_properties_df(gtf_path)
exon_properties = GRanges(exon_properties_df)
names(exon_properties) = exon_properties$exon_row

#include in output, cols in defaults but can be overwritten
#exons_properties_df_cols = c("exon_id","to_coordinate","gene_id","gene_name","gene_type","exon_number","exon_row")

## Harmonize
seqlevelsStyle(svs_gr)="UCSC"
seqlevelsStyle(exon_properties)="UCSC"

exon_overlaps= get_reciprocal_overlap_pairs(svs_gr,exon_properties,reciprocal_overlap = 0,svtype_matching = FALSE)

exon_overlaps = exon_overlaps %>% 
  left_join(exon_properties_df[,exons_properties_df_cols],by=c("set2"="exon_row"))  %>%
  left_join(svs_df[,svs_df_overlap_cols],by=c("set1"="patient_sv_name")) %>% unique() 


write.table(exon_overlaps,sv_exon_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)
