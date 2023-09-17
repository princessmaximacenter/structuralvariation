## Fusion sq
## Set flags based on output from script that performs the
## Pairwise overlap with Population SV databases
## Also add start/end of genes and exons
## 
## Last update 2021-12-05
## External databases
## Use svs union output


# Config is loaded upfront 

if(FALSE){
  #source(paste0("~/PycharmProjects/structuralvariation/sv_functional_analysis/default.conf"))
  
  #using config or paths
  
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.conf")
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.docker.conf")
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.conf")
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/run/wilms_v2_20210923/wilms_v2_20210923.PMCID772AAJ.conf")
  
  #alternatively specify 
  #svs_path=""
  #patient_id=""
  
  sv_analysis_type="somatic"
  
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/utils/annotate_pairwise_overlaps.R")
  
  source('/hpc/pmc_gen/ivanbelzen/case_studies/all_germline/PMRBM000AHT_PMRBM000AHV_PMRBM000AHU.conf');
  source('/hpc/pmc_gen/ivanbelzen/case_studies/all_germline/trio_analysis.conf');
  
}

#default
if(!exists("sv_analysis_type")){
  sv_analysis_type="union"
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

## local override todo refactor
make_sv_bp_gene_partnered = function(sv_bp_gene_overlaps, gene_overlap_cols=c("gene_type","gene_name","gene_id"),svs_df_cols = c("patient_sv_name","patient_sv_merged","partner_sv_name","svtype")) {
  
  svs_unique_cols=c(svs_df_cols,gene_overlap_cols)
  svs_unique_cols=svs_unique_cols[svs_unique_cols %in% names(sv_bp_gene_overlaps)]
  
  svs_to_genes = sv_bp_gene_overlaps[,svs_unique_cols] %>% unique()
  
  
  ## For CTX only => rest will not have partner sv name matching like this, 
  # but do not use     filter(patient_sv_name %in% svs_to_genes$partner_sv_name) because then no longer svs with 1 edge 
  #if you dont remove sv bp orientation you will get inflated counts
  
  sv_bp_gene_partnered_ctx = svs_to_genes %>% filter(svtype=="CTX") %>%
    left_join(svs_to_genes %>% select(-svtype),
              by=c("partner_sv_name"="patient_sv_name","patient_sv_name"="partner_sv_name"))
  
  sv_bp_gene_partnered_ctx = sv_bp_gene_partnered_ctx %>% 
    dplyr::rename_with(.cols = ends_with(".y"), function(x){paste0("partner_",substr(x,0,stop = (str_length(x)-2)))}) %>% 
    dplyr::rename_with(.cols = ends_with(".x"), function(x){substr(x,0,stop = (str_length(x)-2))})
  
  sv_bp_gene_partnered_ctx = sv_bp_gene_partnered_ctx %>% dplyr::rename(partner_sv_merged = partner_patient_sv_merged)
  
  #with sv bp orientation for intra-svs
  svs_to_genes_orientation = sv_bp_gene_overlaps[,c(svs_unique_cols,"sv_breakpoint_orientation")] %>% unique()
  
  #in case only 'end' coordinate, include it in the first datafame and will not have a partner in the 2nd
  sv_bp_gene_partnered_intra = svs_to_genes_orientation %>% filter(svtype!="CTX") %>%
    filter(sv_breakpoint_orientation=="start" | !patient_sv_name %in% filter(svs_to_genes_orientation,sv_breakpoint_orientation=="start")$patient_sv_name ) %>% 
    select(-sv_breakpoint_orientation) %>% unique() 
  sv_bp_gene_partnered_intra_2 = svs_to_genes_orientation %>% filter(svtype!="CTX") %>%
    filter(sv_breakpoint_orientation=="end" & patient_sv_name %in% filter(svs_to_genes_orientation,sv_breakpoint_orientation=="start")$patient_sv_name ) %>% 
    select(patient_sv_name,all_of(gene_overlap_cols)) %>% unique() %>%
    dplyr::rename_with(.cols = all_of(gene_overlap_cols), function(x){paste0("partner_",x)}) 
  
  sv_bp_gene_partnered_intra = sv_bp_gene_partnered_intra %>% left_join(sv_bp_gene_partnered_intra_2)
  
  sv_bp_gene_partnered_intra$partner_sv_merged = sv_bp_gene_partnered_intra$patient_sv_merged
  
  sv_bp_gene_partnered_cols = c(names(sv_bp_gene_partnered_ctx),names(sv_bp_gene_partnered_intra)) %>% unique()
  
  sv_bp_gene_partnered = rbind(sv_bp_gene_partnered_ctx[,sv_bp_gene_partnered_cols],
                               sv_bp_gene_partnered_intra[,sv_bp_gene_partnered_cols])
  
  
  #check if all rows are there
  #nrow(sv_bp_gene_overlaps %>% filter(!patient_sv_name %in% sv_bp_gene_partnered$patient_sv_name))!=0
  
  return(sv_bp_gene_partnered)
} 

## Fill path templates
map_template_vars=c('${resources_dir}'=resources_dir,'${merged_svs_dir}'=merged_svs_dir,'${utils_output_dir}'=utils_output_dir,
                    '${cohort_wdir}'=cohort_wdir,'${cohort_identifier}'=cohort_identifier,'${patient_basename}'=patient$basename)

#input
if(grepl("somatic",sv_analysis_type)) {
  svs_union_path_template = svs_union_somatic_path_template
}

svs_union_path = stri_replace_all_fixed(svs_union_path_template,names(map_template_vars), map_template_vars,vectorize=F)

sv_bp_gene_overlaps_path = stri_replace_all_fixed(sv_bp_gene_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)
sv_bp_exon_overlaps_path = stri_replace_all_fixed(sv_bp_exon_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)
sv_bp_repeatmasker_overlaps_path = stri_replace_all_fixed(sv_bp_repeatmasker_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)
sv_bp_segmental_duplications_overlaps_path = stri_replace_all_fixed(sv_bp_segmental_duplications_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)
#the population overlaps are in the function load_sv_population_db_overlaps()

## output
if(grepl("somatic",sv_analysis_type)) {
  svs_union_anno_path_template = svs_union_anno_somatic_path_template 
  svs_union_anno_multitool_path_template =  svs_union_anno_somatic_multitool_path_template
  
  svs_union_gene_overlaps_anno_path_template = svs_union_gene_overlaps_anno_somatic_path_template 
  svs_union_gene_overlaps_anno_path_multitool_template = svs_union_gene_overlaps_anno_somatic_path_multitool_template 
  
}

svs_union_anno_path = stri_replace_all_fixed(svs_union_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F)
svs_union_anno_multitool_path = stri_replace_all_fixed(svs_union_anno_multitool_path_template,names(map_template_vars), map_template_vars,vectorize=F)

## gene centric
svs_union_gene_overlaps_anno_path = stri_replace_all_fixed(svs_union_gene_overlaps_anno_path_template,names(map_template_vars), map_template_vars,vectorize=F)
svs_union_gene_overlaps_anno_multitool_path = stri_replace_all_fixed(svs_union_gene_overlaps_anno_path_multitool_template,names(map_template_vars), map_template_vars,vectorize=F)



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

if(is.null(svs_df$patient_id) & !is.null(patient_id)){
  svs_df$patient_id=patient_id
  
  #temp code to fix bug with missing $patient_id 
  svs_df$patient_sv_merged = paste0(svs_df$patient_id,"_",svs_df$sv_merged)
  svs_df$partner_sv_name = paste0(svs_df$patient_id,"_",svs_df$partner)
  ##
}
if(is.null(svs_df$patient_sv_merged) & !is.null(patient_id)){
  svs_df$patient_sv_merged = paste0(svs_df$patient_id,"_",svs_df$sv_merged)
}

if(is.null(svs_df$partner_sv_name) & !is.null(patient_id)){
  svs_df$partner_sv_name = paste0(svs_df$patient_id,"_",svs_df$partner)
}

svs_gr = get_svs_gr(svs_df)

#note,  multi tool support is over EVERYTHING right now
#for the sv burden it was after filtering gene overlaps /population svs
multi_tool_support = get_multi_tool_support(svs_df)



### Annotate with population SVS 

sv_databases_lst = c("nstd166","nstd186","dgv") #ignore clinvar
database_overlaps = load_sv_population_db_overlaps(map_template_vars, sv_databases_lst)

#subset overlaps to loaded svs
database_overlaps= database_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
database_overlaps = filter_population_sv_db_overlaps(database_overlaps)


for( sv_database_identifier in sv_databases_lst) {
  overlaps_colname=paste0("flag_sv_db_",sv_database_identifier)
  svs_df=flag_sv_overlap(svs_df,filter(database_overlaps,sv_database==sv_database_identifier),
                         sv_key_col="patient_sv_name",
                         overlaps_key_col="patient_sv_name",overlaps_colname=overlaps_colname)
  
  
  overlaps_colname=paste0("sv_db_",sv_database_identifier,"_overlaps")
  svs_df = annotate_sv_overlap_attr(svs_df,
                           filter(database_overlaps,sv_database==sv_database_identifier),
                           attr_col = "sv_db_bp_name",
                           sv_key_col="patient_sv_name",
                           overlaps_key_col="patient_sv_name",
                           overlaps_colname=overlaps_colname)
}

svs_df = svs_df %>% mutate(flag_sv_population = (flag_sv_db_nstd166 | flag_sv_db_nstd186 | flag_sv_db_dgv))

#checks head(svs_df %>% filter(flag_sv_population))


## Annotate with repeats start/end
#note for the CTX start/end are not providing info on the partner
if(length(Sys.glob(sv_bp_repeatmasker_overlaps_path))==1) { 
  #Load sv repeat overlaps 
  # available cols: repeatmasker_df_cols
  sv_bp_repeatmasker_overlaps = read.table(sv_bp_repeatmasker_overlaps_path,header=T,sep="\t")
  sv_bp_repeatmasker_overlaps = sv_bp_repeatmasker_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, repeat_id=set2) 
  
  #subset to loaded svs
  sv_bp_repeatmasker_overlaps = sv_bp_repeatmasker_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
  
  sv_bp_repeatmasker_overlaps = sv_bp_repeatmasker_overlaps %>% 
    left_join(svs_df[,c("patient_sv_name","patient_sv_merged","partner_sv_name")])#,by="patient_sv_name")
  
  #use partnering instead of a loop with start/end orientation 
  ## this allows you to annotate CTX on partner just like the the genes
  ## flag if any overlap with repeat
  ## flag if both start/end overlap same repeat family
  
  sv_bp_repeatmasker_partnered = make_sv_bp_gene_partnered(sv_bp_repeatmasker_overlaps,gene_overlap_cols = c("repClass"))
  
  sv_bp_repeatmasker_wide = sv_bp_repeatmasker_partnered %>% 
    group_by(patient_sv_name) %>% summarize(start_repeat_class = toString(unique(sort(repClass))),
                                            end_repeat_class = toString(unique(sort(partner_repClass)))) 
  
  sv_bp_repeatmasker_wide = sv_bp_repeatmasker_wide %>% as.data.frame() %>% unique()
  
  svs_df = svs_df %>% 
    left_join(sv_bp_repeatmasker_wide)
  
  col_name_base="flag_repeat"
  svs_df=flag_sv_overlap(svs_df,sv_bp_repeatmasker_overlaps,
                           sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                           overlaps_colname=col_name_base)

  #both breakpoints same repeat class 
  svs_df=flag_sv_overlap(svs_df,
                         sv_bp_repeatmasker_partnered %>% filter(repClass==partner_repClass),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=paste0(c(col_name_base,"both_bp"),collapse = "_"))
  
} else {
  print(paste0("WARNING: missing sv bp repeatmasker overlaps for patient ",patient$patient_id," : ",sv_bp_repeatmasker_overlaps_path))
  #quit()
}

if(length(Sys.glob(sv_bp_segmental_duplications_overlaps_path))==1) { 
  #Load sv segdup overlaps 
  # available cols: segmental_duplications_df_cols
  sv_bp_segmental_duplications_overlaps = read.table(sv_bp_segmental_duplications_overlaps_path,header=T,sep="\t")
  sv_bp_segmental_duplications_overlaps = sv_bp_segmental_duplications_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, segdup_id=set2) 
  
  #subset to loaded svs
  sv_bp_segmental_duplications_overlaps = sv_bp_segmental_duplications_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
  
  sv_bp_segmental_duplications_overlaps = sv_bp_segmental_duplications_overlaps %>% 
    left_join(svs_df[,c("patient_sv_name","patient_sv_merged","partner_sv_name")])#,by="patient_sv_name")
  
  sv_bp_segmental_duplications_partnered = make_sv_bp_gene_partnered(sv_bp_segmental_duplications_overlaps,gene_overlap_cols = c("sd_name"))
  
  if(FALSE) {
    #not easily processed
    sv_bp_segmental_duplications_wide = sv_bp_segmental_duplications_partnered %>% 
    group_by(patient_sv_name) %>% summarize(start_segdup = toString(unique(sort(sd_name))),
                                            end_segdup = toString(unique(sort(partner_sd_name)))) 
  
    sv_bp_segmental_duplications_wide = sv_bp_segmental_duplications_wide %>% as.data.frame() %>% unique()
    
    svs_df = svs_df %>% 
      left_join(sv_bp_segmental_duplications_wide)
    
  }
  
  
  col_name_base="flag_segdup"
  svs_df=flag_sv_overlap(svs_df,sv_bp_segmental_duplications_overlaps,
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=col_name_base)
  
  svs_df=flag_sv_overlap(svs_df,
                         sv_bp_segmental_duplications_partnered %>% filter(!is.na(sd_name) & !is.na(partner_sd_name)),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname=paste0(c(col_name_base,"both_bp"),collapse = "_"))
  
  
} else {
  print(paste0("WARNING: missing sv bp segmental duplication overlaps for patient ",patient$patient_id," : ",sv_bp_segmental_duplications_overlaps_path))
  #quit()
}



### Annotate with genes start/end => use function for partnered 

if(length(Sys.glob(sv_bp_gene_overlaps_path))==1) { 
  #Load sv gene overlaps 
  sv_bp_gene_overlaps = read.table(sv_bp_gene_overlaps_path,header=T,sep="\t")
  sv_bp_gene_overlaps = sv_bp_gene_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, gene_id=set2) 
  
  #subset to loaded svs
  sv_bp_gene_overlaps = sv_bp_gene_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
  
  sv_bp_gene_overlaps = sv_bp_gene_overlaps %>% 
    left_join(svs_df[,c("patient_sv_name","patient_sv_merged","partner_sv_name")])#,by="patient_sv_name")
  
  sv_bp_gene_partnered = make_sv_bp_gene_partnered(sv_bp_gene_overlaps)
  
  #sv_bp_gene_overlaps %>% filter(patient_sv_name=="PMRBM000AHT_PMRBM000AHV_PMRBM000AHU_DUP00009680_bp1--DUP00009680_bp2") %>% group_by(patient_sv_name,sv_breakpoint_orientation) %>% summarize(gene_names = toString(unique(sort(gene_name)))) %>% pivot_wider(names_from="sv_breakpoint_orientation",values_from="gene_names")%>% dplyr::rename(genes_start = start, genes_end=end)
  #same with advantage that you also get the CTX partners
  sv_bp_gene_wide = sv_bp_gene_partnered %>% 
    group_by(patient_sv_name) %>% summarize(start_gene_name = toString(unique(sort(gene_name))),
                                            end_gene_name = toString(unique(sort(partner_gene_name)))) 
  
  sv_bp_gene_wide = sv_bp_gene_wide %>% as.data.frame() %>% unique()
  
  svs_df = svs_df %>% 
    left_join(sv_bp_gene_wide)
  
  #Flags for easy filtering
  
  svs_df=flag_sv_overlap(svs_df,sv_bp_gene_overlaps,
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_bp_in_gene")
  
  svs_df=flag_sv_overlap(svs_df,filter(sv_bp_gene_overlaps,gene_type=="protein_coding"),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_bp_in_gene_coding")
  
} else {
 print(paste0("WARNING: missing sv bp gene overlaps for patient ",patient$patient_id," : ",sv_bp_gene_overlaps_path))
  #quit()
}

if(length(Sys.glob(sv_bp_exon_overlaps_path))==1) { 
  
  #Load sv exon overlaps 
  sv_bp_exon_overlaps = read.table(sv_bp_exon_overlaps_path,header=T,sep="\t")
  sv_bp_exon_overlaps = sv_bp_exon_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype,exon_row=set2) 
  
  sv_bp_exon_overlaps = sv_bp_exon_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
  
  sv_bp_exon_overlaps = sv_bp_exon_overlaps %>% 
    left_join(svs_df[,c("patient_sv_name","patient_sv_merged","partner_sv_name")])#,by="patient_sv_name")
  
  sv_bp_exon_partnered = make_sv_bp_gene_partnered(sv_bp_exon_overlaps)
  
  sv_bp_exon_wide = sv_bp_exon_partnered %>% 
    group_by(patient_sv_name) %>% summarize(start_exon_name = toString(unique(sort(gene_name))),
                                            end_exon_name = toString(unique(sort(partner_gene_name)))) 
  
  
  sv_bp_exon_wide = sv_bp_exon_wide %>% as.data.frame() %>% unique()
  
  svs_df = svs_df %>% 
    left_join(sv_bp_exon_wide)
  
  
  svs_df=flag_sv_overlap(svs_df,sv_bp_exon_overlaps,
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_bp_in_exon")
  
  svs_df=flag_sv_overlap(svs_df,filter(sv_bp_exon_overlaps,gene_type=="protein_coding"),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_bp_in_exon_coding")
  
} else {
    print(paste0("WARNING: missing sv bp exon overlaps for patient ",patient$patient_id," : ",sv_bp_exon_overlaps_path))
    #quit()
}


## add gene overlaps to max of 5

sv_gene_overlaps_path= stri_replace_all_fixed(sv_gene_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)

if(length(Sys.glob(sv_gene_overlaps_path))==1) { 
  
  sv_gene_overlaps = read.table(sv_gene_overlaps_path,header=T,sep="\t")
  sv_gene_overlaps = sv_gene_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, gene_id=set2) 
  sv_gene_overlaps = sv_gene_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
 
  svs_df=flag_sv_overlap(svs_df,sv_gene_overlaps,
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_overlap_gene")
  
  svs_df=flag_sv_overlap(svs_df,filter(sv_gene_overlaps,gene_type=="protein_coding"),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_overlap_gene_coding")
  
  
  svs_df = annotate_sv_overlap_attr(svs_df,
                                    sv_gene_overlaps,
                                    attr_max_cnt = 5,
                                    report_overlap_frac = F,
                                    report_attr_cnt=T,
                                    attr_col = "gene_name",
                                    sv_key_col="patient_sv_name",
                                    overlaps_key_col="patient_sv_name",
                                    overlaps_colname="gene_names")
  svs_df[is.na(svs_df$gene_names_cnt),c("gene_names_cnt")]=0
  
  svs_df = annotate_sv_overlap_attr(svs_df,
                                    sv_gene_overlaps,
                                    attr_max_cnt = 5,
                                    report_overlap_frac = T,
                                    report_attr_cnt=F,
                                    attr_col = "gene_name",
                                    sv_key_col="patient_sv_name",
                                    overlaps_key_col="patient_sv_name",
                                    overlaps_colname="gene_overlaps")
  
  
 } else { 
  print(paste0("WARNING: missing sv gene overlaps for patient ",patient$patient_id," : ",sv_gene_overlaps_path))
  #quit()
}

sv_exon_overlaps_path= stri_replace_all_fixed(sv_exon_overlaps_path_template,names(map_template_vars), map_template_vars,vectorize=F)

if(length(Sys.glob(sv_exon_overlaps_path))==1) { 
  
  sv_exon_overlaps = read.table(sv_exon_overlaps_path,header=T,sep="\t")
  sv_exon_overlaps = sv_exon_overlaps %>% dplyr::rename(patient_sv_name = set1, svtype=set1_svtype, exon_row=set2) 
  sv_exon_overlaps = sv_exon_overlaps %>% select(patient_sv_name,gene_name,gene_id,gene_type) %>% unique()
  sv_exon_overlaps = sv_exon_overlaps %>% filter(patient_sv_name %in% svs_df$patient_sv_name)
  
  svs_df=flag_sv_overlap(svs_df,sv_exon_overlaps,
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_overlap_exon")
  
  svs_df=flag_sv_overlap(svs_df,filter(sv_exon_overlaps,gene_type=="protein_coding"),
                         sv_key_col ="patient_sv_name",overlaps_key_col="patient_sv_name",
                         overlaps_colname="flag_overlap_exon_coding")

  svs_df = annotate_sv_overlap_attr(svs_df,
                                    sv_exon_overlaps,
                                    attr_max_cnt = 5,
                                    report_overlap_frac = F,
                                    report_attr_cnt=T,
                                    attr_col = "gene_name",
                                    sv_key_col="patient_sv_name",
                                    overlaps_key_col="patient_sv_name",
                                    overlaps_colname="gene_exon_names")
  svs_df[is.na(svs_df$gene_exon_names_cnt),c("gene_exon_names_cnt")]=0
  
} else { 
  print(paste0("WARNING: missing sv exon overlaps for patient ",patient$patient_id," : ",sv_exon_overlaps_path))
  #quit()
}


###
## TODO elsewhere add logic like 'intronic' / if it affects exon
## add cancer genes? or is that too version specific
## svs likely affecting exons
#svs_df = svs_df %>% mutate(flag_affects_exon = (svtype=="CTX" & flag_overlap_gene==T) | flag_overlap_exon==T )

#get supporting tool count
svs_df = svs_df %>% left_join(multi_tool_support)
svs_df[is.na(svs_df$tool_cnt),c("tool_cnt")]=1

write.table(svs_df,svs_union_anno_path,quote = FALSE,sep = "\t",row.names=FALSE)


#high confidence only / multitool

write.table(svs_df %>% filter(tool_cnt>1),svs_union_anno_multitool_path,quote = FALSE,sep = "\t",row.names=FALSE)



## Gene centric 

#already has svs_df_overlap_cols
#additional to interpret 

#rename overlap cols to be more clear
## see above patient_sv_name = set1, "from"


if(length(Sys.glob(sv_exon_overlaps_path))==0 | 
   length(Sys.glob(sv_gene_overlaps_path))==0 | 
   length(Sys.glob(sv_bp_gene_overlaps_path))==0 | 
   length(Sys.glob(sv_bp_exon_overlaps_path))==0) { 
  print(paste0("EXITING: Skipping gene centric analysis; one of the gene/exon overlap files missing for patient ",patient$patient_id))
  print(sv_exon_overlaps_path)
  print(sv_gene_overlaps_path)
  print(sv_bp_gene_overlaps_path)
  print(sv_bp_exon_overlaps_path)
  quit()
}

 if(nrow(sv_gene_overlaps)==0) {
   print(paste0("EXITING: Skipping gene centric analysis; 0 sv gene overlaps for patient ",patient$patient_id))
   quit()
 }       

sv_gene_overlaps=sv_gene_overlaps %>% 
  dplyr::rename(overlap_sv_gene_frac = overlap_set1_set2, overlap_gene_sv_frac = overlap_set2_set1,
    sv_coordinate = from_coordinate, gene_coordinate = to_coordinate)

#remove cols 
sv_gene_overlaps = sv_gene_overlaps %>% select(-from,-to,-bp_names_origin,-sv_name)

svs_df_anno_cols = c("variant_type","patient_sv_merged","sv_merged_coordinate","tool_cnt",
  "flag_sv_population","sv_db_nstd166_overlaps","sv_db_nstd186_overlaps","sv_db_dgv_overlaps",
  "start_repeat_class","end_repeat_class","flag_repeat","flag_repeat_both_bp","flag_segdup","flag_segdup_both_bp")

svs_df_anno_cols = svs_df_anno_cols[svs_df_anno_cols %in% names(svs_df)]
sv_gene_overlaps = sv_gene_overlaps %>% left_join(svs_df[,c("patient_sv_name",svs_df_anno_cols)])

#add flags specific to gene-sv pairs 
sv_gene_overlaps$gene_sv_identifier = paste0(sv_gene_overlaps$patient_sv_name,"_",sv_gene_overlaps$gene_id)

#exon overlaps
sv_exon_overlaps$gene_sv_identifier = paste0(sv_exon_overlaps$patient_sv_name,"_",sv_exon_overlaps$gene_id)

sv_gene_overlaps=flag_sv_overlap(sv_gene_overlaps,sv_exon_overlaps,
                       sv_key_col ="gene_sv_identifier",overlaps_key_col="gene_sv_identifier",
                       overlaps_colname="flag_overlap_exon")

#sv bp gene
sv_bp_gene_overlaps$gene_sv_identifier = paste0(sv_bp_gene_overlaps$patient_sv_name,"_",sv_bp_gene_overlaps$gene_id)

sv_gene_overlaps=flag_sv_overlap(sv_gene_overlaps,sv_bp_gene_overlaps,
                                 sv_key_col ="gene_sv_identifier",overlaps_key_col="gene_sv_identifier",
                                 overlaps_colname="flag_bp_in_gene")

#sv bp exon
sv_bp_exon_overlaps$gene_sv_identifier = paste0(sv_bp_exon_overlaps$patient_sv_name,"_",sv_bp_exon_overlaps$gene_id)

sv_gene_overlaps=flag_sv_overlap(sv_gene_overlaps,sv_bp_exon_overlaps,
                                 sv_key_col ="gene_sv_identifier",overlaps_key_col="gene_sv_identifier",
                                 overlaps_colname="flag_bp_in_exon")


#multitool support specific to that gene - merged sv pair 
gene_overlap_multi_tool_support = sv_gene_overlaps  %>% 
  group_by(patient_sv_merged,gene_id) %>% summarize(gene_overlap_tool_cnt=length(unique(tool)))

## NOTE: this does not work 100% as intended, because if you have 2 svs and 1 overlaps with exon, than the one with no overlap will still get a '1' if you left join afterwards. 
#maybe just resolve it by forcing all that are false for the flag as 0 to prevent confusion
#or just accept that the tally is a property of the sv merged not of the sv name. therefore if the sv itself is false for the flag but still has 1 then something is going on.

exon_overlap_multi_tool_support = sv_exon_overlaps  %>% 
  left_join(svs_df[,c("patient_sv_name","patient_sv_merged","tool")]) %>%
  group_by(patient_sv_merged,gene_id) %>% summarize(exon_overlap_tool_cnt=length(unique(tool)))

gene_sv_bp_overlap_multi_tool_support = sv_bp_gene_overlaps  %>% 
  left_join(svs_df[,c("patient_sv_name","patient_sv_merged","tool")]) %>%
  group_by(patient_sv_merged,gene_id) %>% summarize(gene_sv_bp_overlap_tool_cnt=length(unique(tool)))
 
exon_sv_bp_overlap_multi_tool_support = sv_bp_exon_overlaps  %>% 
  left_join(svs_df[,c("patient_sv_name","patient_sv_merged","tool")]) %>%
  group_by(patient_sv_merged,gene_id) %>% summarize(exon_sv_bp_overlap_tool_cnt=length(unique(tool)))

sv_gene_overlaps = sv_gene_overlaps %>% 
  left_join(gene_overlap_multi_tool_support,by=c("patient_sv_merged","gene_id")) %>% 
  left_join(exon_overlap_multi_tool_support,by=c("patient_sv_merged","gene_id")) %>%
  left_join(gene_sv_bp_overlap_multi_tool_support,by=c("patient_sv_merged","gene_id")) %>% 
  left_join(exon_sv_bp_overlap_multi_tool_support,by=c("patient_sv_merged","gene_id"))

#sv_gene_overlaps[is.na(sv_gene_overlaps$gene_overlap_tool_cnt), c("gene_overlap_tool_cnt")]=0
#sv_gene_overlaps[is.na(sv_gene_overlaps$exon_overlap_tool_cnt), c("exon_overlap_tool_cnt")]=0
#sv_gene_overlaps[is.na(sv_gene_overlaps$gene_sv_bp_overlap_tool_cnt), c("gene_sv_bp_overlap_tool_cnt")]=0
#sv_gene_overlaps[is.na(sv_gene_overlaps$exon_sv_bp_overlap_tool_cnt), c("exon_sv_bp_overlap_tool_cnt")]=0


write.table(sv_gene_overlaps,svs_union_gene_overlaps_anno_path,quote = FALSE,sep = "\t",row.names=FALSE)


#high confidence only / multitool
write.table(sv_gene_overlaps %>% filter(tool_cnt>1),svs_union_gene_overlaps_anno_multitool_path,quote = FALSE,sep = "\t",row.names=FALSE)


