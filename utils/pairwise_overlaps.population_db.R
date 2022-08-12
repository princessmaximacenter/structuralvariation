## Pairwise overlap 

# Config is loaded upfront 

## Update 2022-02-15 use svs ranges instead of union, as intermediate step after processing before merging
if(FALSE){
  #using config or paths
  
  source("~/structuralvariation/sv_functional_analysis/default.conf")
  source("~/structuralvariation/sv_functional_analysis/default.docker.conf")
  source("~/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf")
  source("~/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.patient.conf")
#  source("~/structuralvariation/utils/pairwise_overlaps.population_db.R")
  
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
source(paste0(utils_script_dir,"functions.population_svs.R"))

## Fill path templates
map_template_vars=c('${resources_dir}'=resources_dir,'${merged_svs_dir}'=merged_svs_dir,'${utils_output_dir}'=utils_output_dir,
                    '${cohort_wdir}'=cohort_wdir,'${cohort_identifier}'=cohort_identifier,'${patient_basename}'=patient$basename)

svs_ranges_path = stri_replace_all_fixed(svs_ranges_path_template,names(map_template_vars), map_template_vars,vectorize=F)

if(!exists("svs_path")){
  svs_path=svs_ranges_path
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


#for comparison with databases
seqlevelsStyle(svs_gr)="Ensembl"

## For the CTX and other single breakpoints: resize if <30 
svs_gr[width(svs_gr)<30] = flank(svs_gr[width(svs_gr)<30],width = 15,both=T)


#####
## Per datdabase find and save overlaps 

## NCBI multiple databases with standarized formats
databases_lst = c("nstd102","nstd166","nstd186")#,"estd219")

#nstd102 ClinVar: Clinical Structural Variants)  
#nstd186 (NCBI Curated Common Structural Variants)
#nstd166 (gnomAD Structural Variants)
#estd219 (1000 Genomes Consortium Phase 3 Integrated SV) note: very large

#sv_database_identifier="nstd166"
for( sv_database_identifier in databases_lst) {

  db_map_template_vars = c(map_template_vars,'${sv_database_identifier}'=sv_database_identifier)
  sv_database_path = stri_replace_all_fixed(sv_database_path_template,names(db_map_template_vars), db_map_template_vars,vectorize=F)
  database_overlaps_path = stri_replace_all_fixed(sv_population_db_overlaps_path_template,names(db_map_template_vars), db_map_template_vars,vectorize=F)
  
  if(length(Sys.glob(sv_database_path))!=1) {
    print(paste0("WARNING: sv database not found ",sv_database_identifier, " (",sv_database_path,") skipping..."))
    next()
  }
  sv_database_gr = load_sv_database_vcf(sv_database_path)
  
  ##prepare sv database
  mcols(sv_database_gr)$svtype = sv_database_gr$SVTYPE
  sv_database_ranges = make_range_sv_database(sv_database_gr) 
  
  ## cant do SV type matching during overlaps because not exactly the same
  database_overlaps = get_reciprocal_overlap_pairs(sv_database_ranges,svs_gr,reciprocal_overlap = 0.5,svtype_matching = FALSE)
  
  if(nrow(database_overlaps)==0) {
    next()
  }
  
  database_overlaps$sv_database=sv_database_identifier
  
  ## some annotation to make interpretation easier
  
  sv_database_df = as.data.frame(sv_database_ranges[unique(database_overlaps$set1)])
  sv_database_df$to_coordinate = paste0("chr",sv_database_df$seqnames,":",sv_database_df$start,"-",sv_database_df$end)
  sv_database_df = sv_database_df %>% dplyr::rename(sv_db_svlen = SVLEN, sv_db_af = AF, sv_db_ac=AC )
  
  database_overlaps = database_overlaps %>%  
    left_join(sv_database_df[,sv_database_df_cols],by=c("set1"="bp_name"))  %>%
    left_join(svs_df[,svs_df_overlap_cols],by=c("set2"="patient_sv_name")) %>% unique() 
  
  
  write.table(database_overlaps,database_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)
  
}


## DGV database

## Another type of datbase
sv_database_identifier ="dgv"
db_map_template_vars = c(map_template_vars,'${sv_database_identifier}'=sv_database_identifier)

dgv_db_path = stri_replace_all_fixed(dgv_db_path_template,names(db_map_template_vars), db_map_template_vars,vectorize=F)
if(length(Sys.glob(dgv_db_path))!=1) {
  print(paste0("WARNING: sv database not found ",sv_database_identifier, " (",dgv_db_path,") skipping..."))
  quit()
}


database_overlaps_path = stri_replace_all_fixed(sv_population_db_overlaps_path_template,names(db_map_template_vars), db_map_template_vars,vectorize=F)

dgv_db = read.table(dgv_db_path,header=T,sep="\t",stringsAsFactors=F)
dgv_db = dgv_db %>% dplyr::rename(bp_name = variantaccession, sv_db_reference = reference, seqnames=chr, svtype=varianttype, sv_db_variantsubtype=variantsubtype)

dgv_db$svtype = as.character(dgv_db$svtype)
dgv_db[grepl("duplication",dgv_db$variantsubtype),c("svtype")]="DUP"
dgv_db[dgv_db$variantsubtype %in% c("gain"),c("svtype")]="DUP"

dgv_db[grepl("deletion",dgv_db$variantsubtype),c("svtype")]="DEL"
dgv_db[dgv_db$variantsubtype %in% c("loss"),c("svtype")]="DEL"

dgv_db[grepl("insertion",dgv_db$variantsubtype),c("svtype")]="INS"
dgv_db[grepl("inversion",dgv_db$variantsubtype),c("svtype")]="INV"

dgv_db[dgv_db$variantsubtype %in% c("gain+loss","complex"),c("svtype")]="other"

dgv_gr = GRanges(dgv_db[dgv_db$seqnames != "",])
names(dgv_gr) = dgv_gr$bp_name

database_overlaps = get_reciprocal_overlap_pairs(dgv_gr,svs_gr,reciprocal_overlap = 0.5,svtype_matching = FALSE)

if(nrow(database_overlaps)==0) {
  quit()
}
database_overlaps$sv_database=sv_database_identifier

dgv_db = dgv_db[dgv_db$bp_name %in% unique(database_overlaps$set1),]
dgv_db$to_coordinate = paste0("chr",dgv_db$seqnames,":",dgv_db$start,"-",dgv_db$end)

database_overlaps = database_overlaps %>%  
  left_join(dgv_db[,dgv_db_cols],by=c("set1"="bp_name")) %>%
  left_join(svs_df[,svs_df_overlap_cols],by=c("set2"="patient_sv_name")) %>% unique() 


write.table(database_overlaps,database_overlaps_path,quote = F,sep = "\t",row.names=F,col.names = T)





