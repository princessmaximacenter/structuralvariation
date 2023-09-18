## make circos plots
## Last update: 2022-03-04
library(circlize)
library(tidyverse, quietly=TRUE)
library(stringr, quietly=TRUE)
library(stringi)
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

source(paste0(wdir,"solids.conf"))


## Outputs

## Inputs

#to select segments to display
log2foldchange_cutoff=0.2 
maf_cutoff=0.4
autosomes= c(paste("chr",1:22,sep=""),"chrX")


intra_svtype_col= data.frame("INV"="purple","DEL"="blue","DUP"="red")  %>% t() %>% as.data.frame() 
intra_svtype_col$svtype=rownames(intra_svtype_col)
colnames(intra_svtype_col)=c("color","svtype")

cna_call_cols=c("seqnames","start","end","call")
cna_l2fc_cols=c("seqnames","start","end","cr_l2fc_50")
cna_maf_cols=c("seqnames","start","end","maf_50")

results_dir=paste0(cohort_results_dir,"circos/")



## Read in cohort data ----



cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]

patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

automated_baseline_shifts = read.table(automated_baseline_shifts_path,sep="\t",header=T)
manual_baseline_shifts = read.table(manual_baseline_shifts_path,sep="\t",header=T)
baseline_shifts= get_baseline_correction(cohort,automated_baseline_shifts,manual_baseline_shifts)

cancer_type_abbrev = read.table(cancer_type_abbrev_path,sep="\t",header=T)
cohort = cohort %>% left_join(cancer_type_abbrev)

#RESOURCES----

map_template_vars=c('${resources_dir}'=resources_dir)
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband



## TODO refactor per patient based on pipeline output
#load svs 
#filter svs
#merge svs
#also needs cytoband/partner chrom



## Load low tumor AF SVs ---

no_taf_svs_df = load_cohort_svs(cohort,c('${merged_svs_dir}'=merged_svs_dir),svs_path_template=svs_union_anno_somatic_multitool_path_template,flag_add_patient_label = T)

no_taf_svs_df = no_taf_svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
no_taf_svs_df = no_taf_svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
no_taf_svs_df = annotate_sv_multitool_support(no_taf_svs_df,sv_tumor_af_threshold = 0)

no_taf_merged_svs=make_merged_svs(no_taf_svs_df)
no_taf_merged_svs = no_taf_merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  

no_taf_merged_svs = no_taf_merged_svs %>% filter(flag_in_filtered_svs) %>% filter(patient_label %in% cohort$patient_label)
no_taf_merged_svs = no_taf_merged_svs %>% filter(chrom %in% autosomes)
resolve_orphan_svs(no_taf_merged_svs,no_taf_svs_df)

no_taf_merged_svs_gr = GRanges(no_taf_merged_svs$sv_merged_coordinate)
mcols(no_taf_merged_svs_gr) = no_taf_merged_svs 
names(no_taf_merged_svs_gr) = no_taf_merged_svs_gr$patient_sv_merged


## Load SVs ----


unfiltered_svs_df = load_cohort_svs(cohort,c('${merged_svs_dir}'=merged_svs_dir),svs_path_template=svs_union_anno_somatic_multitool_path_template,flag_add_patient_label = T)
unfiltered_svs_df = unfiltered_svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
unfiltered_svs_df = unfiltered_svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
unfiltered_svs_df = annotate_sv_multitool_support(unfiltered_svs_df,sv_tumor_af_threshold = 0.1)
unfiltered_svs_df = unfiltered_svs_df %>% unique()
#if(length(Sys.glob(merged_svs_path)) ==1 & flag_load_from_file==T) {
#merged_svs = read.table(merged_svs_path,sep="\t", header=T)
#} else {

merged_svs = make_merged_svs(unfiltered_svs_df)
merged_svs = merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  
#export merged svs? > in omics 
#}

merged_svs = merged_svs %>% filter(flag_in_filtered_svs) %>% filter(patient_label %in% cohort$patient_label) %>% unique()
#sets globals therefore no output
resolve_orphan_svs(merged_svs,unfiltered_svs_df)

#check consistency merged svs and svs df 
merged_svs %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% nrow() == 0
unfiltered_svs_df %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% filter(patient_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() == 0






## add complex svs
merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.20230701.tsv")
merged_svs_classes = read.table(merged_svs_classes_path,sep="\t",header=T)

merged_svs = merged_svs %>% left_join(merged_svs_classes %>% select(patient_sv_merged,contains("complex")))

## Amplicons ----
amplicon_merged_segments_path = paste0(cohort_results_dir,"amplicon_merged_segments.20230701.tsv")
amplicons_merged_seg_df = read.table(amplicon_merged_segments_path,sep = "\t",header=T)
amplicons_merged_seg = GRanges(amplicons_merged_seg_df$cna_coordinate)
mcols(amplicons_merged_seg) = amplicons_merged_seg_df[,c("patient_label","patient_cna_merged_id","cr_l2fc_50","footprint_id")]
names(amplicons_merged_seg) = amplicons_merged_seg$patient_cna_merged_id


no_taf_map_amplicons_merged_svs = get_reciprocal_overlap_pairs(resize_gr_distance(amplicons_merged_seg,1e3),no_taf_merged_svs_gr,reciprocal_overlap = 0,svtype_matching = F)
no_taf_map_amplicons_merged_svs = no_taf_map_amplicons_merged_svs %>% dplyr::rename(patient_cna_merged_id=set1,patient_sv_merged=set2) 

no_taf_map_amplicons_merged_svs = no_taf_map_amplicons_merged_svs  %>% select(-from,-to) %>% 
  left_join(amplicons_merged_seg_df,by="patient_cna_merged_id") %>% 
  left_join(no_taf_merged_svs[,c("patient_label","patient_sv_merged","sv_merged_coordinate")],by="patient_sv_merged") %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y)


no_taf_map_amplicons_merged_svs_start_end = get_reciprocal_overlap_pairs_start_end(no_taf_merged_svs_gr,resize_gr_distance(amplicons_merged_seg,1e3),reciprocal_overlap = 0,svtype_matching = F)
no_taf_map_amplicons_merged_svs_start_end = no_taf_map_amplicons_merged_svs_start_end %>% dplyr::rename(patient_sv_merged=set1,patient_cna_merged_id=set2) 

no_taf_map_amplicons_merged_svs_start_end = no_taf_map_amplicons_merged_svs_start_end  %>% select(-from,-to) %>% 
  left_join(amplicons_merged_seg_df,by="patient_cna_merged_id") %>% 
  left_join(no_taf_merged_svs[,c("patient_label","patient_sv_merged","sv_merged_coordinate")],by="patient_sv_merged") %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y)


no_taf_map_amplicons_merged_svs_start_end = no_taf_map_amplicons_merged_svs_start_end %>% select(c("patient_label","patient_sv_merged","patient_cna_merged_id","sv_breakpoint_orientation")) %>% pivot_wider(id_cols = c("patient_label","patient_sv_merged"),names_from="sv_breakpoint_orientation",values_from = "patient_cna_merged_id")
no_taf_map_amplicons_merged_svs_start_end = no_taf_map_amplicons_merged_svs_start_end %>% mutate(both_sides=!is.na(start)&!is.na(end))

## Run circos ----

flag_scale_amp=F
flag_no_taf=F
flag_amplicons_ctx_both_sides=T #both sides of ctx need to be in amplicon, disable to getinsertion points
plot_type="full_genome"
plot_type="amplicons"
plot_type="complex_svs"
plot_type="reciprocal_ctx"
plot_type="complex_svs_amplicons"
flag_color_by_complex_sv=T

for(pid in cohort$patient_label){
  patient = filter(cohort,patient_label==pid)
  print(paste0("Running: ",patient$patient_label))

  
  segments_path = paste0(cna_data_dir,patient$tumor_id,"*_WGS*",cna_seg_file_ext)
  segments_file=Sys.glob(segments_path)
  if(length(segments_file)!=1){ 
    print("WARNING: input unclear")
    print(segments_path)
    next()
  }

  modeled_seg = read_modeled_seg(segments_file)
  modeled_seg$tumor_id = patient$tumor_id
  modeled_seg = modeled_seg %>% left_join(baseline_shifts[,c("tumor_id","baseline_correction")])
  modeled_seg[is.na(modeled_seg$baseline_correction),c("baseline_correction")]=0
  modeled_seg$cr_l2fc_50 = modeled_seg$cr_l2fc_50+modeled_seg$baseline_correction
  modeled_seg = modeled_seg %>% call_cna(cna_cutoff=log2foldchange_cutoff, maf_cutoff=maf_cutoff)
  modeled_seg = modeled_seg %>% filter(seqnames %in% autosomes)

  plot_seg = modeled_seg[modeled_seg$call!="0",cna_call_cols]

  scaled_seg=modeled_seg[,cna_l2fc_cols]
  if(flag_scale_amp) { 
    scaled_seg[scaled_seg$cr_l2fc_50>0.6,c("cr_l2fc_50")]=0.6 #for unscaled gain comment out
  }
  scaled_seg[scaled_seg$cr_l2fc_50<(-0.9),c("cr_l2fc_50")]=(-0.9)
  
  maf_seg=modeled_seg[!is.na(modeled_seg$maf_50),cna_maf_cols]
  maf_seg$maf_50=0.5-maf_seg$maf_50
  
  #svs
  if(flag_no_taf) {
    patients_svs_df = no_taf_merged_svs %>% filter(patient_label==pid)
  } else {
    patients_svs_df = merged_svs %>% filter(patient_label==pid)
  }
  
  if(plot_type=="amplicons") {
    patients_svs_df = patients_svs_df %>% 
      filter(patient_label==pid & patient_sv_merged %in% filter(no_taf_map_amplicons_merged_svs_start_end,both_sides)$patient_sv_merged)
    if(flag_amplicons_ctx_both_sides) {
      patients_svs_df = patients_svs_df %>% 
      filter(svtype!="CTX" |  partner_sv_merged %in% filter(no_taf_map_amplicons_merged_svs_start_end,both_sides)$patient_sv_merged)
    }

  }
  
  
  if(plot_type=="complex_svs") {
    complex_sv_lst=unique(filter(patients_svs_df,flag_is_complex&!is.na(complex_sv_id))$complex_sv_id)
    complex_sv_lst=unique(filter(patients_svs_df,!is.na(complex_sv_id))$complex_sv_id)
    if(length(complex_sv_lst) ==0 ) { 
      print("Patient has no complex svs, skipping") 
      next()
    }
  }
  
  if(plot_type=="complex_svs_amplicons") {
    #all amplicon type complex svs
    
    if(flag_no_taf) {
    amplicon_svs_df = patients_svs_df %>% 
      filter(patient_label==pid & 
               patient_sv_merged %in% filter(no_taf_map_amplicons_merged_svs_start_end,both_sides)$patient_sv_merged)
    if(flag_amplicons_ctx_both_sides) {
      amplicon_svs_df = amplicon_svs_df %>% 
        filter(svtype!="CTX" |  
                 partner_sv_merged %in% filter(no_taf_map_amplicons_merged_svs_start_end,both_sides)$patient_sv_merged)
    }
    } else {
      amplicon_svs_df=data.frame()
    }
    complex_svs_df = merged_svs_classes %>% filter(patient_label==pid & complex_sv_class_super=="amplicon")
    
    patients_svs_df = rbind_no_colmatch(amplicon_svs_df,complex_svs_df)
    
  }
    
  if(plot_type=="reciprocal_ctx") {
  complex_sv_lst=unique(filter(patients_svs_df,complex_sv_class=="reciprocal_ctx")$complex_sv_id)
  if(length(complex_sv_lst) ==0 ) { 
    print("Patient has no reciprocal_ctx, skipping") 
    next()
  }
  }
  display_sv_cols=c("patient_sv_merged","sv_merged_coordinate","svtype")
  if(plot_type=="complex_svs" | plot_type=="reciprocal_ctx" | plot_type=="complex_svs_amplicons") {
    display_sv_cols = c("complex_sv_id",display_sv_cols)
  }
  display_ctx = patients_svs_df %>% filter(svtype=="CTX") %>%
    select(all_of(c(display_sv_cols,"partner_sv_merged","partner_sv_merged_coordinate"))) %>%
    filter(!is.na(sv_merged_coordinate)&!is.na(partner_sv_merged_coordinate)) %>% unique()
  
  intra_svs = patients_svs_df %>% filter(svtype!="CTX") %>% select(all_of(display_sv_cols))
  
  display_intra_svs = GRanges(intra_svs$sv_merged_coordinate) %>% as.data.frame() 
  display_intra_svs$patient_sv_merged=intra_svs$patient_sv_merged
  display_intra_svs = display_intra_svs %>% left_join(intra_svs) %>% left_join(intra_svtype_col,by="svtype")
  
  if(flag_color_by_complex_sv) {
    complex_sv_colors = rand_color(length(unique(patients_svs_df$complex_sv_id)), transparency = 0.5,luminosity = "bright")
    complex_sv_colors_df = cbind(unique(patients_svs_df$complex_sv_id),complex_sv_colors) %>% as.data.frame()
    colnames(complex_sv_colors_df) = c("complex_sv_id","color")
    display_intra_svs = display_intra_svs %>% select(-color) %>% left_join(complex_sv_colors_df)
    display_ctx = display_ctx  %>% left_join(complex_sv_colors_df)
    
  }
  
  
  pdf_path = paste0(results_dir,patient$patient_label,".",patient$cancer_type_abbrev,".svs")
  if(flag_no_taf==F){
    pdf_path = paste0(pdf_path,".taf01")
  }
  if(flag_scale_amp==F) {
    pdf_path = paste0(pdf_path,".unscaled_gain")
  }
  
  
  ## Full genome
  if(plot_type=="full_genome") {
    pdf(paste0(pdf_path,".pdf"))
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
    circos.clear()
    
    
    init_circos_default()
    draw_cna_call(plot_seg,draw_loh=T)
    draw_ctx(display_ctx)
    draw_intra_svs(display_intra_svs)
    circos.clear()
    
    dev.off()
  
  } else if(plot_type=="amplicons" | plot_type=="complex_svs_amplicons") {
      #pdf
    amplicon_transparancy = 0.9
    
    if(flag_amplicons_ctx_both_sides==F) {
      pdf_path = paste0(pdf_path,".ctx_single_side")
    }
    if(plot_type=="complex_svs_amplicons") {
      pdf_path = paste0(pdf_path,".complex_svs")
      amplicon_transparancy = 0.85
    } 
    pdf(paste0(pdf_path,".amplicons.pdf"))
    
    patients_svs_df_selection=patients_svs_df
    #chr_lst = unique(patients_svs_df_selection$chrom)
    chr_lst = unique(c(patients_svs_df_selection$chrom,filter(patients_svs_df_selection,svtype=="CTX")$partner_chrom))
    chr_lst=chr_lst[gtools::mixedorder(chr_lst)]
    #scaled_seg[scaled_seg$cr_l2fc_50>2,c("cr_l2fc_50")]=2 #moderate scaling
    
    init_circos_default(chr_lst,start_degree=180)
    selected_seg = scaled_seg %>% call_cna() %>% filter(seqnames %in% chr_lst) #separate line bc chrX can be missing
    if(nrow(selected_seg[selected_seg$call!="0",])>0){
      draw_cna_l2fc(selected_seg)
    }
    #draw_cna_maf(maf_seg)
    draw_ctx(display_ctx %>% filter(patient_sv_merged %in% patients_svs_df_selection$patient_sv_merged))
    draw_intra_svs(display_intra_svs  %>% filter(patient_sv_merged %in% patients_svs_df_selection$patient_sv_merged),transparency=amplicon_transparancy)
    
    circos.clear()
    
    dev.off()

   } else if(plot_type=="complex_svs" | plot_type=="reciprocal_ctx") {
    for(cid in complex_sv_lst) {
      pdf(paste0(pdf_path,".",cid,".pdf"))
      #pdf(paste0(pdf_path,".chr17.pdf"))
      #pdf(paste0(pdf_path,".chr2_amplicons.pdf"))
      #pdf(paste0(pdf_path,".chr2_chr4.pdf"))
      
      patients_svs_df_selection = patients_svs_df %>% filter(complex_sv_id==cid)
      #patients_svs_df_selection = patients_svs_df %>% filter(chrom %in% c("chr2") & !is.na(complex_sv_id))#,"chr2","chr13")) #%>% filter(chrom %in% c("chr1","chr2"))

      #patients_svs_df_selection = patients_svs_df %>% filter(chrom %in% c("chr2","chr13"))# & partner_chrom %in% c("chr2","chr13"))
      
      #chr_lst = unique(patients_svs_df_selection$chrom)
      chr_lst = unique(c(patients_svs_df_selection$chrom,filter(patients_svs_df_selection,svtype=="CTX")$partner_chrom))
      chr_lst=chr_lst[gtools::mixedorder(chr_lst)]
      #scaled_seg[scaled_seg$cr_l2fc_50>2,c("cr_l2fc_50")]=2 #moderate scaling
      
      init_circos_default(chr_lst,start_degree=180)
      selected_seg = scaled_seg %>% call_cna() %>% filter(seqnames %in% chr_lst) #separate line bc chrX can be missing
      if(nrow(selected_seg[selected_seg$call!="0",])>0){
      draw_cna_l2fc(selected_seg)
      }
      #draw_cna_maf(maf_seg)
      
      #draw_ctx(display_ctx %>% filter(patient_sv_merged %in% patients_svs_df_selection$patient_sv_merged))
      #draw_intra_svs(display_intra_svs  %>% filter(patient_sv_merged %in% patients_svs_df_selection$patient_sv_merged),transparency=0.8)
      
      draw_ctx(display_ctx %>% filter(complex_sv_id==cid))
      draw_intra_svs(display_intra_svs %>% filter(complex_sv_id==cid),transparency=0.9)
      
      circos.clear()
      
      dev.off()
    } 
  }

}

