## Merge SVs
## per patient, merge svs from tools to find confident ones
# Tools: manta, delly, gridss

## Last update 2022-02-16
## Note: previously merge_svs.somatic.docker.R
### spit up the non-somatic merge pipeline (first parse_sv_vcf.R, filter population overlaps, then merge_svs.R)
## this is somatic only and needs to be run after the population overlaps have been generated. 

#load default.conf 
#optional loading of default.docker.conf 
#load patient conf 
#load unit test conf 

# Config order is important because of overwrites
# general unit test config goes last, which uses patient config for pathnames, 
# also unit test identifier sometimes used in path names 

## Change log 
## Update 2022-02-16
## Refactor, previously merge_svs.somatic.docker.R
## Changed filtering from tumor AF >0.05; to read support and population db overlap
## spit up the non-somatic merge pipeline (first parse_sv_vcf.R, filter population overlaps, then merge_svs.R)
## this is somatic only and needs to be run after the population overlaps have been generated. 

## update 2021-07-05
# Greatly speed up SVs to ranges with vectorisation
# Intermediate output SVs as ranges
# Single function file: functions.svs.R no longer requires functions from fusion-sq, and no longer uses functions.annotate_sv_genes.R 
# Prefilter VCF to remove low qual, and apply threshold on tumor AF. 

## Update 2021-07-08 
# template paths with variable replacement to enable docker overrides more easily without changing a lot of configs
###

## Update 2021-10-07
# somatic only script to see if the merged svs would be different 

## TEST CODE
if(FALSE){

  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.conf');
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.docker.conf');
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.PMCID057AAK.conf');
  source('/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf');
  
}


if(is.null(patient) | is.null(cohort_identifier)) {
  print("EXIT: patient undefined, need config file with tumor/normal ids ")
  quit()
} 


read_support_threshold = 6

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

svs_ranges_somatic_path = stri_replace_all_fixed(svs_ranges_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)
svs_union_somatic_path = stri_replace_all_fixed(svs_union_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)
svs_merged_somatic_path = stri_replace_all_fixed(svs_merged_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)

svs_filtering_log_somatic_path = stri_replace_all_fixed(svs_filtering_log_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)

vcf_manta_somatic_path = stri_replace_all_fixed(vcf_manta_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)
vcf_manta_germline_path = stri_replace_all_fixed(vcf_manta_germline_path_template,names(map_template_vars), map_template_vars,vectorize=F)

vcf_gridss_somatic_path = stri_replace_all_fixed(vcf_gridss_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)
vcf_delly_somatic_path = stri_replace_all_fixed(vcf_delly_somatic_path_template,names(map_template_vars), map_template_vars,vectorize=F)


manta_present=delly_present=gridss_present=T
if(length(Sys.glob(vcf_manta_somatic_path))!=1) {
  print(paste0("Manta somatic missing: ", vcf_manta_somatic_path) )
  manta_present=F
}
if(length(Sys.glob(vcf_manta_germline_path))!=1) {
  print(paste0("Manta diploid missing", vcf_manta_germline_path))
  manta_present=F
}
if(length(Sys.glob(vcf_delly_somatic_path))!=1) {
  print(paste0("Delly somatic missing", vcf_delly_somatic_path))
  delly_present=F
}
if(length(Sys.glob(vcf_gridss_somatic_path))!=1) {
  print(paste0("GRIDSS somatic missing", vcf_gridss_somatic_path))
  gridss_present=F
}


##



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
if(delly_present) vcf_delly_somatic_path = prefilter_pass(vcf_delly_somatic_path)
if(gridss_present) vcf_gridss_somatic_path = prefilter_pass(vcf_gridss_somatic_path)

manta_gr=delly_gr=gridss_gr=GRanges()
if(manta_present) manta_gr = read_manta_sv_vcf(vcf_manta_germline_path, vcf_manta_somatic_path,patient)
if(delly_present) delly_gr = read_delly_sv_vcf(vcf_delly_somatic_path,patient)
if(gridss_present) gridss_gr = read_gridss_sv_vcf(vcf_gridss_somatic_path,patient)


all_gr=c(manta_gr[manta_gr$somatic==TRUE],delly_gr,gridss_gr)
all_gr$somatic=T
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

all_gr = all_gr[all_gr$partner %in% names(all_gr)] #remove unpartnered

## Make partnered ranges
read_support_colnames = colnames(metadata)[grepl("DV|DR|RV|RR|REF|SR|PR|VF",colnames(metadata))]
read_support_colnames  = read_support_colnames[grepl("tumor|normal",read_support_colnames)]      

sv_metadata_cols =  c("sourceId",  "svtype", "svLen", "partner","FILTER","QUAL","ALT","REF",
                      "insLen", "HOMLEN", "tumor_af", "normal_af",  "somatic", "tool", "supporting_reads_tumor", read_support_colnames)

svs = make_range_svs(all_gr,sv_metadata_cols)
svs$patient_sv_name = paste0(patient$patient_id,"_",svs$sv_name)

svs_df = as.data.frame(svs)
svs_df = unique(svs_df)
svs_df$coordinate = paste0(svs_df$seqnames,":",svs_df$start,"-",svs_df$end,":",svs_df$strand)
write.table(svs_df,svs_ranges_somatic_path,sep="\t",col.names = T,row.names = F,quote=F)


cnt_svs_after_support_filter=nrow(svs_df)

#Remove population db overlaps

sv_databases_lst = c("nstd166","nstd186","dgv") #ignore clinvar
database_overlaps = load_sv_population_db_overlaps(map_template_vars, sv_databases_lst)
#subset overlaps to loaded svs to speed up filtering
database_overlaps= database_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
database_overlaps = database_overlaps %>% filter(overlap_set1_set2>0.9 & overlap_set2_set1 > 0.9)
database_overlaps = filter_population_sv_db_overlaps(database_overlaps)

svs_pop_db_df = svs_df %>% filter(patient_sv_name %in% database_overlaps$patient_sv_name )


svs_df = svs_df %>% filter(!patient_sv_name %in% database_overlaps$patient_sv_name )
svs = svs[svs$patient_sv_name %in% svs_df$patient_sv_name]

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
overlap_merged = find_same_sv(svs_ranges,
                              svs_ranges,
                              reciprocal_overlap = 0.5,svtype_matching = T,ignore_strand = F)

## Build supporting SV dataframe
if(length(overlap_merged)>0){
  #hash the sv merged identifiers for uq names
  map_merged_to_sv_names = overlap_merged %>% group_by(sv_merged) %>% 
    summarize(sv_names = toString(unique(sort(c(set1)))),
              sv_merged_hash= paste0("merged_",md5(sv_names)))
  overlap_merged = overlap_merged %>% left_join(map_merged_to_sv_names[,c("sv_merged","sv_merged_hash")],by="sv_merged")
  
  overlap_merged = overlap_merged %>% select(-sv_merged) %>% dplyr::rename(sv_merged = sv_merged_hash)
  
  ## check if sv only assigned once 
  if(nrow(unique(overlap_merged[,c("set1","sv_merged")]))!=length(unique(overlap_merged$set1))) {
    uq_merged = unique(overlap_merged[,c("set1","sv_merged")])
    print( overlap_merged[overlap_merged$set1 %in% uq_merged[duplicated(uq_merged$set1),c("set1")],])
    
    print("WARNING SV assigned multiple times")
    print(patient$patient_identifier)
    break
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

svs_df = annotate_variant_af_class(svs_df)

partnered_coord = svs_df[,c("sv_name","partner","coordinate")] %>%
  left_join(svs_df[,c("sv_name","partner","coordinate")],
            by=c("partner"="sv_name","sv_name"="partner")) %>% 
  dplyr::rename(partner_coordinate= coordinate.y, coordinate=coordinate.x)

svs_df=svs_df %>% left_join(partnered_coord[,c("sv_name","partner","partner_coordinate")],by=c("sv_name","partner"))

svs_df$patient_sv_merged = paste0(svs_df$patient_id,"_",svs_df$sv_merged)
write.table(svs_df,svs_union_somatic_path,sep="\t",col.names = T,row.names = F,quote=F)



## NOTE< almost all code above is borrowed from combine_wgs_support.R

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


write.table(sv_merged_df,svs_merged_somatic_path,sep="\t",col.names = T,row.names = F,quote=F)

multi_tool_support = get_multi_tool_support(svs_df)

svs_no_multitool = svs_df %>% filter(!patient_sv_merged %in% multi_tool_support$patient_sv_merged)
svs_multitool = svs_df %>% filter(patient_sv_merged %in% multi_tool_support$patient_sv_merged)

## Write log for filtering

con <- file(svs_filtering_log_somatic_path, open="wt")
writeLines(paste0("gr_delly: ",length(delly_gr)), con)
writeLines(paste0("gr_manta: ",length(manta_gr)), con)
writeLines(paste0("gr_gridss: ",length(gridss_gr)), con)
writeLines(paste0("read_support_threshold: ",read_support_threshold), con)
writeLines(paste0("sv_bp_insufficient_reads: ",nrow(sv_bp_insufficient_reads)), con)
writeLines(paste0("svs_after_support_filter: ",cnt_svs_after_support_filter), con)
writeLines(paste0("svs_pop_db: ",nrow(svs_pop_db_df)), con)
writeLines(paste0("svs_filtered: ",nrow(svs_df)), con)
writeLines(paste0("svs_multitool: ",nrow(svs_multitool)), con)
writeLines(paste0("svs_no_multitool: ",nrow(svs_no_multitool)), con)
writeLines(paste0("svs_multitool_merged: ",length(unique(svs_multitool$patient_sv_merged))), con)
writeLines(paste0("svs_multitool_merged_taf>0.05: ",length(unique(filter(svs_multitool,!is.na(tumor_af)&tumor_af>0.05)$patient_sv_merged))), con)
close(con)
