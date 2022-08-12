## Annotate modelled segments from GATK CNV pipeline
## get CN ratio and alt-AF to determine CNAs and LOH 
#script: wilms.annotate_cna_chromosome.R

## Last update: 2022-03-10
## enable processing chrX if available => functions file
## Prevent spurious LOH calls for missing MAF => functions file
## Update: 2021-11-06
## Refactor for HPC: run for single patient / input file
#should be able to run for every possible seg input 


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

## settings
expected_autosomal_length = 3088269832-156040895-57227415

##

map_template_vars=c('${resources_dir}'=resources_dir,'${utils_output_dir}'=utils_output_dir,'${biomaterial_id}'=run_biomaterial_id)



#output paths ----
chr_arm_cna_overlap_path = stri_replace_all_fixed(chr_arm_cna_overlap_path_template,names(map_template_vars), map_template_vars,vectorize=F)
cytoband_cna_overlap_path = stri_replace_all_fixed(cytoband_cna_overlap_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chr_arm_cna_path = stri_replace_all_fixed(chr_arm_cna_path_template,names(map_template_vars), map_template_vars,vectorize=F)
cytoband_cna_path = stri_replace_all_fixed(cytoband_cna_path_template,names(map_template_vars), map_template_vars,vectorize=F)
fga_per_chrom_path = stri_replace_all_fixed(fga_per_chrom_path_template,names(map_template_vars), map_template_vars,vectorize=F)


### RESOURCES ----
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband
   
#Read CN data ---- 
  
  

  ## Update 2022-03-10 chrX if available
  autosomes = c(paste("chr",1:22,sep=""))
  contig_lengths = data.frame(seqnames=character(), length=character(), stringsAsFactors = F)
  contig_lengths = get_contig_lengths(segments_file,chromosomes = c(autosomes,"chrX"))
  
  if(sum(contig_lengths[contig_lengths$seqnames %in% autosomes,c("length")]) != expected_autosomal_length) {
    print("WARNING: autosomal length different than expected, is it really hg38?")
  }
  
  modeled_seg= read_modeled_seg(segments_file)

  ## Per chromosome FGA
  fga_per_chrom = data.frame()
  for(chrom in contig_lengths$seqnames) {
    segments = modeled_seg[modeled_seg$seqnames==chrom & abs(modeled_seg$cr_l2fc_50) > 0.2,]

    fga_entry = c()
    fga_entry$seqnames = chrom
    fga_entry$cna_bases = sum(segments$width)
    fga_entry$seg_cnt = nrow(modeled_seg[modeled_seg$seqnames==chrom,])
    fga_entry$seg_cna_cnt = nrow(segments)
    fga_entry$fga = fga_entry$cna_bases / contig_lengths[contig_lengths$seqnames==chrom,c("length")]
    fga_entry$max_cna = ifelse(fga_entry$cna_bases>0,max(segments$cr_l2fc_50),0)
    fga_entry$min_cna = ifelse(fga_entry$cna_bases>0,min(segments$cr_l2fc_50),0)
    
    fga_per_chrom = rbind(fga_per_chrom,fga_entry)
    
  }
  
  ## Overall FGA
  fga_entry = c()
  fga_entry$seqnames = "total"
  fga_entry$cna_bases = sum(fga_per_chrom$cna_bases)
  fga_entry$seg_cnt = sum(fga_per_chrom$seg_cnt)
  fga_entry$seg_cna_cnt = sum(fga_per_chrom$seg_cna_cnt)
  fga_entry$fga = fga_entry$cna_bases / sum(contig_lengths$length)
  fga_entry$max_cna = max(fga_per_chrom$max_cna)
  fga_entry$min_cna = min(fga_per_chrom$min_cna)
  
  fga_per_chrom = rbind(fga_per_chrom,fga_entry)

  if(exists("run_patient_id")) { fga_per_chrom$patient_id=run_patient_id }
  
  write.table(fga_per_chrom,fga_per_chrom_path,sep="\t",quote=F,row.names = F, col.names = T)
  

## Get and store reciprocal overlaps ----
  
  seg_gr=GRanges(modeled_seg)
  names(seg_gr)=seg_gr$cna_id

  #default
  #modeled_seg_cols=c("cna_id","width","cr_l2fc_50","maf_50","cr_l2fc_10","maf_10","cr_l2fc_90","maf_90")  

  cytoband_cna_overlaps = get_reciprocal_overlap_pairs(seg_gr,chromosome_bands,reciprocal_overlap = 0,svtype_matching = F)
  cytoband_cna_overlaps = cytoband_cna_overlaps %>% left_join(modeled_seg[,modeled_seg_cols],by=c("set1"="cna_id"))

  
  write.table(cytoband_cna_overlaps,cytoband_cna_overlap_path,sep="\t",quote=F,row.names = F, col.names = T)
  
  chr_arm_cna_overlaps = get_reciprocal_overlap_pairs(seg_gr,chr_arms,reciprocal_overlap = 0,svtype_matching = F)
  chr_arm_cna_overlaps = chr_arm_cna_overlaps %>% left_join(modeled_seg[,modeled_seg_cols],by=c("set1"="cna_id"))
  
  write.table(chr_arm_cna_overlaps,chr_arm_cna_overlap_path,sep="\t",quote=F,row.names = F, col.names = T)
  
  
  
## Summarize over chr arms ----
    chr_arm_cna_overlaps = read.table(chr_arm_cna_overlap_path,sep="\t",header=T)
  
    chr_arm_cna_overlaps = chr_arm_cna_overlaps %>% dplyr::rename(chr_arm=set2, cna_id=set1)
    chr_arm_cna_overlaps = chr_arm_cna_overlaps %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4)
    
    chr_arm_cna = get_overlaps_modeled_seg_summary(chr_arm_cna_overlaps,group_cols=c("chr_arm"))
    ## Note that : if MAF is NA everywhere then becomes 0 and then later on causes LOH call. -> reset afterwards to not mess up the function
    
    chr_arm_cna_extra = chr_arm_cna_overlaps %>% 
      group_by(chr_arm) %>%
      summarize(maf_50_sd = sd(maf_50,na.rm=T),
                cr_l2fc_50_sd = sd(cr_l2fc_50,na.rm=T),
                .groups="drop") %>% as.data.frame()


    chr_arm_cna = chr_arm_cna %>% left_join(chr_arm_cna_extra)
    
  #stability by percentage of bases 
  threshold=0.33 ## cr threshold 1.33 and 0.67, for MAF only max because often jumping around between 0 and 0.5 
  
  thresh = chr_arm_cna %>% 
    dplyr::rename(mean_cr_l2fc_50=cr_l2fc_50,mean_maf_50=maf_50) %>% 
    dplyr::select(chr_arm,mean_cr_l2fc_50,mean_maf_50) %>%
    mutate(min_cr = ifelse(abs(mean_cr_l2fc_50)<0.1, -0.1, (1-threshold)*mean_cr_l2fc_50),
           max_cr = ifelse(abs(mean_cr_l2fc_50)<0.1, 0.1, (1+threshold)*mean_cr_l2fc_50),
          # min_maf = (1-(threshold/2))*mean_maf_50,
           max_maf = max(0.1,(1+(threshold/2))*mean_maf_50))
  
  #min max should be opposite if negative numbers 
  
  chr_arm_cr_stability = chr_arm_cna_overlaps %>% left_join(thresh,by="chr_arm") %>%
    group_by(chr_arm,mean_cr_l2fc_50) %>%
    filter( 
      (min_cr<max_cr & cr_l2fc_50 > min_cr & cr_l2fc_50 < max_cr) | 
      (max_cr<min_cr & cr_l2fc_50 < min_cr & cr_l2fc_50 > max_cr)
    ) %>%
    summarize(seg_covered_cr_stable = sum(overlap_set2_set1))
  
  chr_arm_maf_stability = chr_arm_cna_overlaps %>% left_join(thresh,by="chr_arm") %>%
    group_by(chr_arm,mean_maf_50) %>%
    filter(  maf_50 < max_maf) %>%
    summarize(seg_covered_maf_stable = sum(overlap_set2_set1))
  
  
  chr_arm_cna = chr_arm_cna  %>% 
    left_join(chr_arm_cr_stability[,c("chr_arm","seg_covered_cr_stable")],by="chr_arm") %>% 
    left_join(chr_arm_maf_stability[,c("chr_arm","seg_covered_maf_stable")],by="chr_arm")
  
  ## how to do this for neutral chromosomes with bit of noise? => the min/max 0.1
  ## mean_cr_l2fc_50 is misleading if low seg covered but acen/stalk/gvar excluded anyway
 
  chr_arm_cna = chr_arm_cna %>% as.data.frame() %>% call_cna()
  
  chr_arm_cna = chr_arm_cna %>% call_cn_stability() %>% call_cn_alteration() 
  
  #add patient id but do not rely on it being there?.. 
  if(exists("run_patient_id")) { chr_arm_cna$patient_id=run_patient_id }
  write.table(chr_arm_cna,chr_arm_cna_path,sep="\t",quote=F,row.names = F, col.names = T)
  
  
### Summarize over cytobands ----
  
    cytoband_cna_overlaps = read.table(cytoband_cna_overlap_path,sep="\t",header=T)
    
    cytoband_cna_overlaps = cytoband_cna_overlaps %>% dplyr::rename(cytoband=set2, cna_id=set1)
    cytoband_cna_overlaps = cytoband_cna_overlaps %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4) 
      
    cytoband_cna = get_overlaps_modeled_seg_summary(cytoband_cna_overlaps,group_cols=c("cytoband"))

    #prevent false LOH calls
    cytoband_cna[cytoband_cna$seg_covered_maf_data==0,c("maf_50","maf_90","maf_10")] = NA
    
    cytoband_cna = cytoband_cna %>% as.data.frame() %>% call_cna()
    #check cytoband_cna  %>% filter(grepl("chrX",cytoband))
    
    #add patient id but do not rely on it being there?.. 
    if(exists("run_patient_id")) { cytoband_cna$patient_id=run_patient_id }
    
    write.table(cytoband_cna,cytoband_cna_path,sep="\t",quote=F,row.names = F, col.names = T)
 


