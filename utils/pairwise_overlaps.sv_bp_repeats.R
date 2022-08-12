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
sv_bp_repeatmasker_overlaps_path = stri_replace_all_fixed(sv_bp_repeatmasker_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)
sv_bp_segmental_duplications_overlaps_path = stri_replace_all_fixed(sv_bp_segmental_duplications_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)

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
#                "sv_name","bp_names_origin","patient_id","partner","tool","tumor_af","normal_af")

#as fail safe, only report columns that exist
svs_df_overlap_cols = svs_df_overlap_cols[svs_df_overlap_cols %in% names(svs_df)]

### RESOURCES 
repeatmasker_path = stri_replace_all_fixed(repeatmasker_path_template,names(map_template_vars), map_template_vars,vectorize=F)
segmental_duplications_path = stri_replace_all_fixed(segmental_duplications_path_template,names(map_template_vars), map_template_vars,vectorize=F)

##repeatmasker
repeatmasker = read.table(repeatmasker_path,header=T,sep="\t",comment.char = "")
repeatmasker = repeatmasker %>% dplyr::rename("seqnames"="genoName","start"="genoStart","end"="genoEnd")
repeatmasker = repeatmasker %>% filter(repClass %in% c("LINE","SINE","LTR") & abs(repLeft)<50)
repeatmasker$repeat_id = paste0("rep_",1:nrow(repeatmasker))
repeatmasker = GRanges(repeatmasker)
names(repeatmasker) = repeatmasker$repeat_id

repeatmasker_df = as.data.frame(repeatmasker)
rownames(repeatmasker_df) = repeatmasker_df$repeat_id

repeatmasker_overlaps = get_reciprocal_overlap_pairs_start_end(svs_gr,repeatmasker,reciprocal_overlap = 0,svtype_matching = FALSE)
repeatmasker_overlaps = repeatmasker_overlaps %>%  
  left_join(repeatmasker_df[,repeatmasker_df_cols],by=c("set2"="repeat_id"))  %>%
  left_join(svs_df[,svs_df_overlap_cols],by=c("set1"="patient_sv_name")) %>% unique() 

write.table(repeatmasker_overlaps,sv_bp_repeatmasker_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)

## seg dups
segmental_duplications = read.table(segmental_duplications_path,header=T,sep="\t",comment.char = "")
segmental_duplications = segmental_duplications %>% dplyr::rename("seqnames"="chrom","start"="chromStart","end"="chromEnd","sd_name"="name")
segmental_duplications$segdup_id = paste0("segdup_",1:nrow(segmental_duplications))
segmental_duplications = GRanges(segmental_duplications)
names(segmental_duplications) = segmental_duplications$segdup_id

segmental_duplications_df = as.data.frame(segmental_duplications)

segmental_duplications_overlaps= get_reciprocal_overlap_pairs_start_end(svs_gr,segmental_duplications,reciprocal_overlap = 0,svtype_matching = FALSE)
segmental_duplications_overlaps = segmental_duplications_overlaps %>% 
  left_join(segmental_duplications_df[,segmental_duplications_df_cols],by=c("set2"="segdup_id"))  %>%
  left_join(svs_df[,svs_df_overlap_cols],by=c("set1"="patient_sv_name")) %>% unique() 

write.table(segmental_duplications_overlaps,sv_bp_segmental_duplications_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)
