## make circos plots
## Last update: 2022-03-04
library(circlize)
library(tidyverse, quietly=TRUE)
library(stringr, quietly=TRUE)
library(stringdist, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)

#if("dplyr" %in% (.packages())){
detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
library(plyr)
library(dplyr)

wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(wdir,"functions.cna.R"))
source(paste0(wdir,"functions.general.R"))
source(paste0(wdir,"functions.svs.R"))
source(paste0(wdir,"functions.svs.local.R"))

source(paste0(utils_script_dir,"functions.circos.R"))

source(paste0(wdir,"wilms.conf"))


## Outputs

## Inputs
cna_data_dir = paste0("~/data/cnv/modeled_seg/")
#cna_seg_file_ext = ".tumor.modelFinal.seg"
cna_seg_file_ext = ".modelFinal.seg"
## TODO some are missing tumor => use .modelFinal.seg
## TODO and should use normal Model as comparison

#to select segments to display
log2foldchange_cutoff=0.2 
maf_cutoff=0.4
autosomes= c(paste("chr",1:22,sep=""),"chrX")


## Read in cohort
cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
cohort$patient_identifier = paste0(cohort$patient_id,"_",cohort$rna_id)
cohort$basename = paste0(cohort$tumor_id,"_",cohort$normal_id)
cohort_anno = read.table(cohort_anno_path,sep = "\t", header=T) %>% arrange(patient_id)
cohort = cohort %>% left_join(cohort_anno )
nrow(cohort)

results_dir="/Users/ianthevanbelzen/results/wilms_v3/circos/"

merged_svs_path = "~/results/wilms_v3/merged_sv.somatic.tsv"
merged_svs = read.table(merged_svs_path,sep="\t",header=T)

#Read somatic SVs
svs_df = read.table("~/results/wilms_v2/cohort_sv.somatic.anno.tsv",sep="\t",header=T)
svs_df$patient_sv_merged = paste0(svs_df$patient_id,"_",svs_df$sv_merged)
svs_df$partner_sv_name = paste0(svs_df$patient_id,"_",svs_df$partner) #for multitool ctx
multitool_rs6 = get_multi_tool_support(svs_df)
cohort_somatic_svs_df = svs_df %>% filter(patient_sv_merged %in% multitool_rs6$patient_sv_merged) %>% filter(patient_id %in% cohort$patient_id)

## filtered
cohort_somatic_svs_df = cohort_somatic_svs_df %>% filter(tumor_af>0.1)
cohort_somatic_svs_df = cohort_somatic_svs_df %>% filter(seqnames %in% autosomes)
## for most coding analyses filter by tumor AF >0.1 and then reassess multi tool support 
filtered_svs_multi_tool_support = get_multi_tool_support(cohort_somatic_svs_df)
cohort_somatic_svs_df = cohort_somatic_svs_df %>% filter(patient_sv_merged %in% filtered_svs_multi_tool_support$patient_sv_merged) 

#replacement for merged tumor af mean instead of all >0.1 tumor AF
#cohort_somatic_svs_df = cohort_somatic_svs_df %>% filter(patient_sv_merged %in% filter(merged_svs,tumor_af_mean>0.1)$patient_sv_merged) 


intra_svtype_col= data.frame("INV"="purple","DEL"="blue","DUP"="red")  %>% t() %>% as.data.frame() 
intra_svtype_col$svtype=rownames(intra_svtype_col)
colnames(intra_svtype_col)=c("color","svtype")

cna_call_cols=c("seqnames","start","end","call")
cna_l2fc_cols=c("seqnames","start","end","cr_l2fc_50")
cna_maf_cols=c("seqnames","start","end","maf_50")


for(pid in cohort$patient_id){
#  for(pid in focus_patients){
  patient = filter(cohort,patient_id==pid)
  print(paste0("Running: ",patient$patient_id))

  
  segments_path = paste0(cna_data_dir,patient$tumor_id,"*_WGS*",cna_seg_file_ext)
  segments_file=Sys.glob(segments_path)
  if(length(segments_file)!=1){ 
    print("WARNING: input unclear")
    print(segments_path)
    next()
  }

  modeled_seg = read_modeled_seg(segments_file)
  modeled_seg = modeled_seg %>% call_cna(cna_cutoff=log2foldchange_cutoff, maf_cutoff=maf_cutoff)
  modeled_seg = modeled_seg %>% filter(seqnames %in% autosomes)
  
  plot_seg = modeled_seg[modeled_seg$call!="0",cna_call_cols]

  scaled_seg=modeled_seg[,cna_l2fc_cols]
  scaled_seg[scaled_seg$cr_l2fc_50>0.6,c("cr_l2fc_50")]=0.6 #for unscaled gain comment out
  scaled_seg[scaled_seg$cr_l2fc_50<(-0.9),c("cr_l2fc_50")]=(-0.9)
  
  maf_seg=modeled_seg[!is.na(modeled_seg$maf_50),cna_maf_cols]
  maf_seg$maf_50=0.5-maf_seg$maf_50
  
  #svs
  patients_svs_df = cohort_somatic_svs_df %>% filter(patient_id==pid)  %>% filter(seqnames %in% autosomes)
  ctx = patients_svs_df %>% filter(svtype=="CTX") 
  ctx = ctx %>% get_merged_svs(return_partnered = T) %>% select(sv_merged,sv_merged_coordinate,partner_sv_merged,partner_sv_merged_coordinate) 
  display_ctx = ctx %>% filter(!is.na(sv_merged_coordinate)&!is.na(partner_sv_merged_coordinate))
  
  intra_svs = patients_svs_df %>% filter(svtype!="CTX") 
  
  intra_svs_merged=intra_svs[,c("sv_merged_coordinate","svtype")] %>% unique() 
  intra_svs = GRanges(intra_svs_merged$sv_merged_coordinate) %>% as.data.frame() 
  intra_svs$svtype=intra_svs_merged$svtype
  
  intra_svs = intra_svs %>% left_join(intra_svtype_col,by="svtype")
  
  display_intra_svs = intra_svs 
  
  dev.off()
  #pdf(paste0(results_dir,patient$patient_id,".svs.taf01.pdf"))
  #pdf(paste0(results_dir,patient$patient_id,".svs.taf01.unscaled_gain.pdf"))
  #pdf(paste0(results_dir,patient$patient_id,".svs.taf01.incl_subclonal.pdf"))
  
  #pdf(paste0(results_dir,patient$patient_id,".svs.pdf"))

  ## Full genome

  
  circos.clear()
  
  init_circos_default()
  draw_cna_l2fc(scaled_seg)
  draw_cna_maf(maf_seg)
  draw_ctx(display_ctx)
  draw_intra_svs(display_intra_svs)
  circos.clear()
  
  
  init_circos_default()
  draw_cna_l2fc(scaled_seg)
  draw_ctx(display_ctx)
  draw_intra_svs(display_intra_svs)
  circos.clear()
  
  
  init_circos_default()
  draw_cna_call(plot_seg,draw_loh=T)
  draw_ctx(display_ctx)
  draw_intra_svs(display_intra_svs)
  circos.clear()
  
  
  dev.off()  
  
  
  
  
}
   