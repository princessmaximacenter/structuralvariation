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
## Update 2022-02-08
## Refactor, previously merge_svs.docker.2.R
## Changed filtering from tumor AF >0.05; to read support and population db overlap
### now spit up into 
# - parse_sv_vcf.R which collects vcf from 3 tools and prefilters by supporting reads -> sv ranges 
# sv ranges first retrieve output with population svs 
# - merge_svs.R takes sv ranges and filter out population svs; and merge remaining svs -> sv merged 
# also: sv merged names with md5 hash, and added partner coordinates 

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
  source('~/structuralvariation/sv_functional_analysis/default.conf');
  source('~/structuralvariation/sv_functional_analysis/default.docker.conf');
  source('~/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.patient.conf');
  source('~/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf');
  
}


if(is.null(patient) | is.null(cohort_identifier)) {
  print("EXIT: patient undefined, need config file with tumor/normal ids ")
  quit()
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
  library(openssl)
})

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.svs.variant_type.R"))
source(paste0(utils_script_dir,"functions.population_svs.R"))

## Fill path templates
## replace data paths here instead of in config to enable overwrites for environments eg local/docker 
map_template_vars=c('${input_vcf_data_dir}'=input_vcf_data_dir,'${merged_svs_dir}'=merged_svs_dir,
                    '${recurrent_output_dir}'=recurrent_output_dir,'${utils_output_dir}'=utils_output_dir,
                    '${cohort_wdir}'=cohort_wdir,
                    '${cohort_identifier}'=cohort_identifier,'${patient_basename}'=patient$basename)

svs_ranges_path = stri_replace_all_fixed(svs_ranges_path_template,names(map_template_vars), map_template_vars,vectorize=F)
svs_union_path = stri_replace_all_fixed(svs_union_path_template,names(map_template_vars), map_template_vars,vectorize=F)
svs_merged_path = stri_replace_all_fixed(svs_merged_path_template,names(map_template_vars), map_template_vars,vectorize=F)

svs_filtering_log_path = stri_replace_all_fixed(svs_filtering_log_path_template,names(map_template_vars), map_template_vars,vectorize=F)


##

# load SV ranges 
if(length(Sys.glob(svs_ranges_path))==0) {
  print(paste0("EXIT: not found SV ranges for ",patient$basename," ", svs_ranges_path))
  quit()
}

  svs_df= read.table(svs_ranges_path,sep="\t",header=T)
  #for the pop db overlap filtering
  svs_df$patient_id=patient$patient_id
  svs_df$patient_sv_name = paste0(patient$patient_id,"_",svs_df$sv_name)
  
#prevent crash of merge
## TODO later investigate this 
  svs_df[is.na(svs_df$variant_type),c("variant_type")]="low_af"
  

## Remove population db overlaps
  
  sv_databases_lst = c("nstd166","nstd186","dgv") #ignore clinvar
  database_overlaps = load_sv_population_db_overlaps(map_template_vars, sv_databases_lst)
  
  #subset overlaps to loaded svs
  #database_overlaps= database_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
  database_overlaps = database_overlaps %>% filter(overlap_set1_set2>0.9 & overlap_set2_set1 > 0.9)
  database_overlaps = filter_population_sv_db_overlaps(database_overlaps)
  
  svs_pop_db_df = svs_df %>% filter(patient_sv_name %in% database_overlaps$patient_sv_name )
  
  svs_df = svs_df %>% filter(!patient_sv_name %in% database_overlaps$patient_sv_name )
  
  
  svs=GRanges(svs_df)
  names(svs)=svs$sv_name
  
## Prepare for merging SVs
## For the CTX and other single breakpoints: resize if <30 
## rename to svs_ranges because want to get back to the original ones later
svs_ranges = svs
svs_ranges[width(svs_ranges)<30] = flank(svs_ranges[width(svs_ranges)<30],width = 15,both=T)

#map partnered ranges to bp 
map_sv_range_bp_name= as.data.frame(svs) %>% dplyr::select(bp_name) %>% 
  separate(col=bp_name,into=c("bp_name_head","bp_name_tail"),sep = "--",remove = F,fill = "left") %>% 
  dplyr::select(bp_name,bp_name_head,bp_name_tail) %>% dplyr::rename(sv_name = bp_name) %>%
  gather(key = "bp_name_orient",value="bp_name",-sv_name) %>% dplyr::select(sv_name,bp_name)

map_sv_range_bp_name[is.na(map_sv_range_bp_name$bp_name),c("bp_name")]=map_sv_range_bp_name[is.na(map_sv_range_bp_name$bp_name),c("sv_name")]
map_sv_range_bp_name = unique(map_sv_range_bp_name)


# Find Same SVs
## which of the SVs overlap >50% and make a merged/reduced SV 

#returns pairwise overlaps between SVs and the sv merged they have in common
# also contains sv merged coordinate and the overlap fractions between each with merged and with eachother
## NB: currently not using the merged SV itself in this pipeline, but the identifier to group similar/same SV events
#overlap_merged = find_same_sv(svs_ranges,
#                              svs_ranges,
#                              reciprocal_overlap = 0.5,svtype_matching = T,ignore_strand = F)

overlap_merged=data.frame()
for(vt in c("tumor_specific","low_af","germline")) {
  overlap_merged_vt = find_same_sv(svs_ranges[svs_ranges$variant_type==vt],
                                   svs_ranges[svs_ranges$variant_type==vt],
                                   reciprocal_overlap = 0.5,svtype_matching = T,ignore_strand = F)
  if(nrow(overlap_merged_vt)==0) {next()}
  overlap_merged_vt$variant_type=vt
  overlap_merged=rbind(overlap_merged,overlap_merged_vt)  
}


## Build supporting SV dataframe
if(length(overlap_merged)>0){
  
  #helper for the hashing, otherwise wrong svs are taken together
  overlap_merged$sv_merged = paste0(overlap_merged$sv_merged,"_",overlap_merged$variant_type)
  
  #hash the sv merged identifiers for uq names
  map_merged_to_sv_names = overlap_merged %>% group_by(sv_merged) %>% 
    summarize(sv_names = toString(unique(sort(c(set1)))),
              sv_merged_hash= paste0("merged_",md5(sv_names)))
  overlap_merged = overlap_merged %>% left_join(map_merged_to_sv_names[,c("sv_merged","sv_merged_hash")],by="sv_merged")
  
  overlap_merged = overlap_merged %>% select(-sv_merged) %>% dplyr::rename(sv_merged = sv_merged_hash)
  
  #add variant type back because convenient
  overlap_merged$sv_merged = paste0(overlap_merged$sv_merged,"_",overlap_merged$variant_type)
  
    
  ## check if sv only assigned once 
  if(nrow(unique(overlap_merged[,c("set1","sv_merged")]))!=length(unique(overlap_merged$set1))) {
    uq_merged = unique(overlap_merged[,c("set1","sv_merged")])
    print( overlap_merged[overlap_merged$set1 %in% uq_merged[duplicated(uq_merged$set1),c("set1")],])
    
    print("WARNING SV assigned multiple times")
    print(patient$patient_identifier)
    #break
  }
  
  svs_df = as.data.frame(svs) %>% 
    left_join(overlap_merged[,c("set1","sv_merged","sv_merged_coordinate","overlap_merged_set1","overlap_set1_merged")], by=c("bp_name"="set1")) 
  
} else {
  svs_df = as.data.frame(svs)
  svs_df$sv_merged = NA
  svs_df$overlap_merged_set1 = NA
  svs_df$overlap_set1_merged = NA
  svs_df$sv_merged_coordinate = NA
}

svs_df = unique(svs_df)
svs_df$coordinate = paste0(svs_df$seqnames,":",svs_df$start,"-",svs_df$end,":",svs_df$strand)


## set NA ranges to bp name
svs_df[is.na(svs_df$sv_merged),c("sv_merged")]=
  svs_df[is.na(svs_df$sv_merged),c("bp_name")]
svs_df[is.na(svs_df$sv_merged_coordinate),c("sv_merged_coordinate")]=
  svs_df[is.na(svs_df$sv_merged_coordinate),c("coordinate")]

#add the 'sv name' also explicitly
svs_df$sv_name=svs_df$bp_name

#add partner coordinates
partnered_coord = svs_df[,c("sv_name","partner","coordinate")] %>%
  left_join(svs_df[,c("sv_name","partner","coordinate")],
            by=c("partner"="sv_name","sv_name"="partner")) %>% 
  dplyr::rename(partner_coordinate= coordinate.y, coordinate=coordinate.x)

svs_df=svs_df %>% left_join(partnered_coord[,c("sv_name","partner","partner_coordinate")],by=c("sv_name","partner"))

svs_df$patient_id=patient$patient_id
svs_df$patient_sv_name = paste0(svs_df$patient_id,"_",svs_df$sv_name)
svs_df$patient_sv_merged = paste0(svs_df$patient_id,"_",svs_df$sv_merged)

write.table(svs_df,svs_union_path,sep="\t",col.names = T,row.names = F,quote=F)


#update 2022-02-15 made consistent with somatic merged => coordinateS and partner_coordinates
sv_merged_df = svs_df %>% dplyr::rename(tumor_af_bp = tumor_af, normal_af_bp = normal_af) %>%
  group_by(sv_merged,sv_merged_coordinate) %>% summarize(
    pass_all = all(FILTER=="PASS"),
    FILTER = toString(unique(sort(FILTER))),
    sv_names = toString(unique(sort(as.character(sv_name)))),
    partners = toString(unique(sort(as.character(partner)))),
    coordinates = toString(unique(sort(coordinate))),
    partner_coordinates = toString(unique(sort(partner_coordinate))),
    
    tools=toString(unique(sort(tool))),
    tools_cnt = length(unique(tool)),
    tumor_af=mean(tumor_af_bp,na.rm=T), normal_af=mean(normal_af_bp,na.rm=T),
    tumor_af_spread=(max(tumor_af_bp,na.rm=T)-min(tumor_af,na.rm=T)), normal_af_spread=(max(normal_af_bp,na.rm=T)-min(normal_af_bp,na.rm=T)), 
    
    variant_type=toString(unique(sort(variant_type))),
    
    svtype=toString(unique(sort(svtype))), 
    svlen=ifelse(!grepl("CTX",svtype),mean(svLen,na.rm=T),NA),
    svlen_spread=ifelse(!is.na(svlen),(max(svLen,na.rm=T)-min(svLen,na.rm = T)),NA),
    .groups="keep") %>% ungroup() 


## show partnered bps

partnered = sv_merged_df[,c("sv_merged","partners","sv_names")] %>% 
  left_join(sv_merged_df[,c("sv_merged","partners","sv_names")],
            by=c("partners"="sv_names","sv_names"="partners")) %>% 
  dplyr::rename(partner_sv_merged= sv_merged.y, sv_merged=sv_merged.x)

sv_merged_df = sv_merged_df %>% left_join(partnered[,c("sv_merged","partner_sv_merged")],by = "sv_merged")


write.table(sv_merged_df,svs_merged_path,sep="\t",col.names = T,row.names = F,quote=F)


multi_tool_support = get_multi_tool_support(svs_df)

svs_no_multitool = svs_df %>% filter(!patient_sv_merged %in% multi_tool_support$patient_sv_merged)
svs_multitool = svs_df %>% filter(patient_sv_merged %in% multi_tool_support$patient_sv_merged)

## Write log for filtering

con <- file(svs_filtering_log_path, open="a")
writeLines(paste0("svs_pop_db: ",nrow(svs_pop_db_df)), con)
writeLines(paste0("svs_filtered: ",nrow(svs_df)), con)
writeLines(paste0("svs_multitool: ",nrow(svs_multitool)), con)
writeLines(paste0("svs_no_multitool: ",nrow(svs_no_multitool)), con)
writeLines(paste0("svs_multitool_merged: ",length(unique(svs_multitool$patient_sv_merged))), con)
writeLines(paste0("svs_multitool_merged_taf>0.05: ",length(unique(filter(svs_multitool,!is.na(tumor_af)&tumor_af>0.05)$patient_sv_merged))), con)
close(con)
