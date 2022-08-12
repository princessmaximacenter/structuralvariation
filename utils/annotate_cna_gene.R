## Annotate modelled segments from GATK CNV pipeline
## get CN ratio and alt-AF to determine CNAs and LOH 
### script: annotate_cna_gene.R 

## Last update: 2022-03-10
## enable processing chrX if available

## Note: 2021-21-21
## Refactor for HPC: run for single patient / input file
#should be able to run for every possible seg input 
### similar to annotate_cna_chromosome.R 


if(FALSE) {
  #local
  wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
  source(paste0(wdir,"default.conf"))
}

if(FALSE){
  #config example 
  
  source("~/structuralvariation/sv_functional_analysis/default.conf")
  source("~/structuralvariation/sv_functional_analysis/default.docker.conf")
  
  run_biomaterial_id="PMABM000GFX" #required
  run_patient_id=""   #add patient id but do not rely on it being there?.. not sure if sensible
  segments_path = "" #sys glob based on biomaterial id if not provided
}

#sys glob with biomaterial id if no path provided 
if(!exists("run_biomaterial_id")) {
  print("No biomaterial id provided, exiting")
  quit()
}

#2022-03 after rerun of pipeline all tumor samples have 'tumor' in filename
#add to default config, then can override with hpc config 
#cna_data_dir = paste0("~/data/cnv/modeled_seg/")

## TODO and should use normal Model as comparison
if(!exists("cna_seg_file_ext")) {
  #cna_seg_file_ext = ".tumor.modelFinal.seg"
  cna_seg_file_ext = ".modelFinal.seg"
}
if(!exists("sequencing_strategy")) {
  sequencing_strategy="WGS"
}
if(!exists("segments_path")) {
  segments_path = paste0(cna_data_dir,run_biomaterial_id,"*_",sequencing_strategy,"*",cna_seg_file_ext)
}

segments_file=Sys.glob(segments_path)
if(length(segments_file)!=1){ 
  print("WARNING: input unclear")
  print(segments_path)
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

map_template_vars=c('${resources_dir}'=resources_dir)

gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)
gene_properties_df = get_gene_properties_df(gtf_path)

#include in output, cols in defaults but can be overwritten
#gene_properties_df_cols = c("gene_id","to_coordinate","gene_name","gene_type","gene_width")
gene_properties = GRanges(gene_properties_df)


## Input
## todo add these as templates and parse with stringi
gene_cna_overlap_path = paste0(utils_output_dir,gene_cna_overlap_outfile,run_biomaterial_id,".tsv")
gene_cna_path = paste0(utils_output_dir,gene_cna_outfile,run_biomaterial_id,".tsv")


  modeled_seg= read_modeled_seg(segments_file)
  modeled_seg$from_coordinate = paste0(modeled_seg$seqnames,":",modeled_seg$start,"-",modeled_seg$end)
  
  #default
  #modeled_seg_cols=c("cna_id","width","cr_l2fc_50","maf_50","cr_l2fc_10","maf_10","cr_l2fc_90","maf_90")  
  modeled_seg_cols=c(modeled_seg_cols,"from_coordinate")
  modeled_seg_cols = modeled_seg_cols[modeled_seg_cols %in% names(modeled_seg)]
  
  seg_gr=GRanges(modeled_seg)
  names(seg_gr)=seg_gr$cna_id

  gene_cna_overlaps = get_reciprocal_overlap_pairs(seg_gr,gene_properties,reciprocal_overlap = 0,svtype_matching = F)
  
  gene_cna_overlaps = gene_cna_overlaps %>% left_join(modeled_seg[,modeled_seg_cols],by=c("set1"="cna_id"))
  gene_cna_overlaps = gene_cna_overlaps %>% left_join(gene_properties_df[,gene_properties_df_cols],by=c("set2"="gene_id"))
  
 # write.table(gene_cna_overlaps,gene_cna_overlap_path,sep="\t",quote=F,row.names = F, col.names = T)

  # summarize over the feature
  gene_cna_overlaps = gene_cna_overlaps %>% dplyr::rename(gene_id=set2, cna_id=set1)
  gene_cna_overlaps = gene_cna_overlaps %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4)
  
  group_cols=c("gene_name","gene_id","gene_type")
  gene_cna = get_overlaps_modeled_seg_summary(gene_cna_overlaps,group_cols)
  
  gene_cna_extra = gene_cna_overlaps %>%
    group_by(across(all_of(group_cols))) %>%
    summarize(overlapping_cna=toString(unique(sort(paste0(cna_id," (",round(cr_l2fc_50,2),")")))),
              overlapping_cna_id=toString(unique(sort(cna_id))),.groups="drop") %>% as.data.frame()
  
  gene_cna = gene_cna %>% left_join(gene_cna_extra)
  
  gene_cna = gene_cna %>% call_cna(cna_cutoff = 0.2,maf_cutoff = 0.4,cna_colname = "cr_l2fc_50",maf_colname = "maf_50")
  
  if(exists("run_patient_id")) { gene_cna$patient_id=run_patient_id }
  
  write.table(gene_cna,gene_cna_path,sep="\t",quote=F,row.names = F, col.names = T)
  

