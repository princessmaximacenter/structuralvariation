## Analyse if amplicons connected by SVs
## Last update: 2023-09-05
# input amplicon merged seg and footprints 
# do (low AF) SVs connect amplicons on different chromosomes or multiple peaks within a large footprint.


suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  #library(stringdist, quietly=TRUE)
})


wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(wdir,"functions.cna.R"))
source(paste0(wdir,"functions.general.R"))
source(paste0(wdir,"functions.svs.R"))

source(paste0(wdir,"solids.conf"))

# Settings ----




## Inputs
run_date=20230701
ctx_chrom_cn_properties_path = paste0(cohort_results_dir,"ctx_chrom_cn_properties.",run_date,".tsv")
unbalanced_ctx_windows_path = paste0(cohort_results_dir,"unbalanced_ctx_windows.",run_date,".tsv")
map_complex_sv_cn_change_path = paste0(cohort_results_dir,"map_complex_sv_cn_change.",run_date,".tsv")


## Outputs
override_if_exists=F
chromosome_cn_calls_path = paste0(cohort_results_dir,"chromosome_level_cn_calls.",run_date,".tsv")
overlap_aneuploidy_unb_ctx_path = paste0(cohort_results_dir,"overlap_aneuploidy_unbalanced_ctx.",run_date,".tsv")
cohort_fga_per_chrom_path = paste0(cohort_results_dir,"cohort_fga_per_chrom.",run_date,".tsv")
cohort_fga_per_chrom_per_call_path = paste0(cohort_results_dir,"cohort_fga_per_chrom_per_call.",run_date,".tsv")
cohort_ploidy_path = paste0(cohort_results_dir,"cohort_ploidy.",run_date,".tsv")
cohort_complex_fga_per_chrom_path = paste0(cohort_results_dir,"cohort_fga_per_chrom.complex.",run_date,".tsv")


# Cohort ----
cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

# Resources ----
## Chromosomes ----
map_template_vars=c('${resources_dir}'=resources_dir)
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)

chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)
chr_arms_df$chr_arm_width=chr_arms_df$width

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband

chr_arm_order = chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")]# %>% flatten_chr()
chrom_order = chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")] %>% unique()# %>% flatten_chr()

chr_arms_unsplit = get_chr_arms(chromosome_bands_df,split_giestain=F)
chr_arms_unsplit = GRanges(chr_arms_unsplit)
names(chr_arms_unsplit) = chr_arms_unsplit$chr_arm
chr_arms_unsplit_df=as.data.frame(chr_arms_unsplit)
chr_arms_unsplit_df = chr_arms_unsplit_df %>% dplyr::rename(chr_arm_width = width)

## NB: cannot remove the acen/gvar regions from chromosomes because it spans the centromere
#start is always 1 but depends on giestain used or not so this is more flexible
chromosomes_df = chr_arms_unsplit_df %>% 
  #  filter(!grepl("acen|gvar|stalk|chrY",chr_arm)) %>% 
  group_by(seqnames) %>% summarize(start = min(start),end=max(end),.groups="drop")
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames
chromosomes_df$chrom_width=chromosomes_df$end-chromosomes_df$start+1 #to make equal with chr arm width



# Read  data ----

## CN data ----
automated_baseline_shifts = read.table(automated_baseline_shifts_path,sep="\t",header=T)
manual_baseline_shifts = read.table(manual_baseline_shifts_path,sep="\t",header=T)
baseline_shifts= get_baseline_correction(cohort,automated_baseline_shifts,manual_baseline_shifts)

segments_df = load_cohort_cn_segments(cohort,map_template_vars=map_template_vars,baseline_shifts=baseline_shifts,flag_add_patient_label=T)
segments_df = segments_df %>% get_gr_coordinate("cna_coordinate")
segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id

## Load unbalanced CTX ----

# 2023-05-05 version using windows chrom, chr arm and 5mb around breakpoints
# updated 2023-07-01 with change in CN data

ctx_chrom_cn_properties =read.table(ctx_chrom_cn_properties_path,header=T,sep="\t")
unbalanced_ctx = ctx_chrom_cn_properties %>% filter(chrom_cn_unbalanced)
unbalanced_ctx = unbalanced_ctx %>%  mutate(chr_arm_orientation = ifelse(grepl("p",chr_arm),"p","q"))

unbalanced_ctx_windows_df =read.table(unbalanced_ctx_windows_path,header=T,sep="\t")

## get chrom/arm fractions ----

chrom_cn_fractions = get_chrom_cn_fractions(segments,segments_df,
                                            chromosomes,chromosomes_df,cohort,ranges_id_col="seqnames",ranges_width_col="chrom_width",
                                            relative_cr = F,return_wide = T)
chrom_cn_fractions = chrom_cn_fractions %>% dplyr::mutate(chrom_mean_cr=cr_l2fc_50)  

chr_arm_cn_fractions = get_chrom_cn_fractions(segments,segments_df,chr_arms,chr_arms_df,cohort,relative_cr = F,return_wide = T,
                                              ranges_id_col="chr_arm",ranges_width_col="chr_arm_width")
chr_arm_cn_fractions = chr_arm_cn_fractions %>% left_join(chr_arms_df[,c("chr_arm","seqnames")])  %>% left_join(cohort[,c("patient_label","cancer_type")])
chr_arm_cn_fractions = chr_arm_cn_fractions %>% filter(!grepl("acen|gvar|stalk",chr_arm))

# FGA ----

## get FGA from chr arm fractions ----

#no acen/gvar/stalk

chr_arm_cn_fractions_long = get_chrom_cn_fractions(segments,segments_df,chr_arms,chr_arms_df,cohort,relative_cr = F,return_wide = F,
                                                   ranges_id_col="chr_arm",ranges_width_col="chr_arm_width")
chr_arm_cn_fractions_long = chr_arm_cn_fractions_long %>% left_join(chr_arms_df[,c("chr_arm","seqnames")])  %>% left_join(cohort[,c("patient_label","cancer_type")])
chr_arm_cn_fractions_long = chr_arm_cn_fractions_long %>% filter(!grepl("acen|gvar|stalk",chr_arm))

cohort_fga_per_chrom=data.frame()
for(pid in unique(chr_arm_cn_fractions_long$patient_label)) {
  fga_per_chrom = chr_arm_cn_fractions_long %>% filter(patient_label==pid) %>%
    dplyr::mutate(width=seg_covered) %>% get_fga_per_chrom(contig_lengths = NULL,flag_use_provided_length = T)
  fga_per_chrom$patient_label=pid
  cohort_fga_per_chrom = rbind(cohort_fga_per_chrom,fga_per_chrom)
}


cohort_fga_per_chrom_per_call=data.frame()

for(cn_type in c("gain","loss","loh")) {
  for(pid in unique(chr_arm_cn_fractions_long$patient_label)) {
    chr_data_patient =  chr_arm_cn_fractions_long %>% filter(patient_label==pid) %>%
      dplyr::mutate(width=seg_covered) 
    tmp_contig_lengths = chr_data_patient %>% group_by(seqnames) %>% summarize(length=sum(width, na.rm=T)) %>% as.data.frame()
    fga_per_chrom =  get_fga_per_chrom(filter(chr_data_patient,call==cn_type),contig_lengths = tmp_contig_lengths)
    fga_per_chrom$patient_label=pid
    fga_per_chrom$cn_type=cn_type
    cohort_fga_per_chrom_per_call = rbind(cohort_fga_per_chrom_per_call,fga_per_chrom)
  }
}


## FGA by complex svs ----

#remove acen/gvar/stalk to calculate proper percentage
map_complex_sv_cn_change = read.table(map_complex_sv_cn_change_path,sep="\t",header = T)
windows_complex_sv_cn_change_gr = GRanges(map_complex_sv_cn_change$window_coordinate)
names(windows_complex_sv_cn_change_gr) = map_complex_sv_cn_change$window_id
mcols(windows_complex_sv_cn_change_gr) = map_complex_sv_cn_change


segments_windows_overlaps_complex = get_segments_windows_overlaps(segments,segments_df,windows_complex_sv_cn_change_gr,
                                                                  windows_df=map_complex_sv_cn_change,window_cols = c("window_id","complex_sv_id"))
segments_windows_overlaps_complex = segments_windows_overlaps_complex %>% select(patient_cna_id,complex_sv_id,patient_label) %>% unique()
#does duplicated exist? yes makes sense and is OK
#segments_windows_overlaps_complex[duplicated(segments_windows_overlaps_complex$patient_cna_id),] 

complex_segments_df = segments_df %>% merge(segments_windows_overlaps_complex)
complex_segments = GRanges(complex_segments_df)
names(complex_segments) = complex_segments$patient_cna_id


#count segs only once for FGA
complex_segments_genome_df = segments_df %>% merge(segments_windows_overlaps_complex %>% select(-complex_sv_id) %>% unique())
complex_segments_genome = GRanges(complex_segments_genome_df)
names(complex_segments_genome) = complex_segments_genome$patient_cna_id


complex_chr_arm_cn_fractions_long = get_chrom_cn_fractions(complex_segments_genome,complex_segments_genome_df,
                                                           chr_arms,chr_arms_df,cohort,relative_cr = F,return_wide = F,
                                                           ranges_id_col="chr_arm",ranges_width_col="chr_arm_width",
                                                           cna_id_col="patient_cna_id")
complex_chr_arm_cn_fractions_long = complex_chr_arm_cn_fractions_long %>% left_join(chr_arms_df[,c("chr_arm","seqnames")])  %>% left_join(cohort[,c("patient_label","cancer_type")])
complex_chr_arm_cn_fractions_long = complex_chr_arm_cn_fractions_long %>% filter(!grepl("acen|gvar|stalk",chr_arm))

complex_cohort_fga_per_chrom=data.frame()
for(pid in unique(complex_chr_arm_cn_fractions_long$patient_label)) {
  chr_data_patient =  chr_arm_cn_fractions_long %>% filter(patient_label==pid) %>% dplyr::mutate(width=seg_covered) 
  tmp_contig_lengths = chr_data_patient %>% group_by(seqnames) %>% summarize(length=sum(width, na.rm=T)) %>% as.data.frame()
  fga_per_chrom =  get_fga_per_chrom(filter(complex_chr_arm_cn_fractions_long,patient_label==pid) %>% dplyr::mutate(width=seg_covered),contig_lengths = tmp_contig_lengths)
  fga_per_chrom$patient_label=pid
  complex_cohort_fga_per_chrom = rbind(complex_cohort_fga_per_chrom,fga_per_chrom)
}


complex_cohort_fga_per_chrom = complex_cohort_fga_per_chrom %>% dplyr::rename(fga_complex=fga) %>%
  left_join(cohort_fga_per_chrom[,c("patient_label","fga","seqnames")],by=c("patient_label","seqnames")) %>%
  mutate(fga_complex_relative=fga_complex/fga) 


## Ploidy  ----

  autosome_chr_arms = chr_arms[!grepl("acen|gvar|stalk|chrY|chrX",chr_arms$chr_arm)]
  autosome_chr_arms$chr_arm_length = width(autosome_chr_arms)
  
  cohort_ploidy = chr_arm_cn_fractions %>% left_join(autosome_chr_arms[,c("chr_arm","chr_arm_length")] %>% as.data.frame()) %>% 
    filter(chr_arm %in% autosome_chr_arms$chr_arm) %>%
    group_by(patient_label) %>% summarize(cr_mean_of_chr_arms= mean(cr_l2fc_50,na.rm=T),
                                          cr_mean_weighted = sum(cr_l2fc_50*chr_arm_length)/sum(autosome_chr_arms$chr_arm_length),
                                          ploidy = 2*(2^(cr_mean_weighted)))




# Get numerical CN changes ----

# 2023-05-24
# 1) Get aneuploidy calls by make_chromosome_cn_calls()
# use selected=T to determine chrom/arm
# 2) Get all unb ctx pass, determine chrom/arm 
# 3) For all unb ctx pass, use p/q arm logic to get the terminal segments. => make into GR
# 4) For all aneuploidy calls, take chrom/arm and call state => make into GR
# 5) overlap aneuploidy and unb ctx GRs. Match on patient label and call state. 
# Those that match are not aneuploidy.
# 

chromosome_cn_calls = make_chromosome_cn_calls(segments,segments_df,chromosomes,chromosomes_df,chr_arms_unsplit,chr_arms_unsplit_df,make_relative_cr = F)
chromosome_cn_calls = chromosome_cn_calls %>% left_join(cohort[,patient_tumor_id_cols]) 

unbalanced_ctx_gr = GRanges(unbalanced_ctx_windows_df)
names(unbalanced_ctx_gr)=unbalanced_ctx_gr$window_id

chromosome_cn_calls = chromosome_cn_calls %>% dplyr::mutate(aneuploidy_id = paste0(patient_label,"_",region,"_"),
                                                            aneuploidy_id = ifelse(region=="chr_arm", paste0(aneuploidy_id,chr_arm),paste0(aneuploidy_id,seqnames)))

chromosome_cn_calls_chr_arm_gr =chromosome_cn_calls %>% filter(region=="chr_arm") %>% left_join(chr_arms_unsplit_df %>% select(chr_arm,seqnames,start,end))  %>% GRanges()
chromosome_cn_calls_chrom_gr =chromosome_cn_calls %>% filter(region=="chrom") %>% left_join(chromosomes_df %>% select(seqnames,seqnames,start,end))  %>% GRanges()

chromosome_cn_calls_gr = c(chromosome_cn_calls_chr_arm_gr,chromosome_cn_calls_chrom_gr)
names(chromosome_cn_calls_gr)=chromosome_cn_calls_gr$aneuploidy_id

aneuploidy_cols = c("aneuploidy_id","patient_label","call","region","selected","cancer_type")
unbalanced_ctx_cols = c("window_id","patient_label","call","region")
overlap_aneuploidy_unb_ctx = get_reciprocal_overlap_pairs(chromosome_cn_calls_gr,unbalanced_ctx_gr,reciprocal_overlap = 0,svtype_matching = F)
overlap_aneuploidy_unb_ctx = overlap_aneuploidy_unb_ctx %>% dplyr::rename(aneuploidy_id=set1,window_id=set2) %>% 
  left_join(chromosome_cn_calls[,aneuploidy_cols],by="aneuploidy_id") %>% 
  left_join(unbalanced_ctx_windows_df[,unbalanced_ctx_cols],by="window_id")

#region does not have to match! I used opposite logic for unbalanced, chosing chromosome over arm because it is not a judgement over the segment but to show the chrom itself is stable.

overlap_aneuploidy_unb_ctx = overlap_aneuploidy_unb_ctx %>% filter(patient_label.x == patient_label.y  & call.x==call.y)
overlap_aneuploidy_unb_ctx = overlap_aneuploidy_unb_ctx %>% filter(overlap_set1_set2>0.7)

chromosome_cn_calls = chromosome_cn_calls  %>% dplyr::mutate(flag_overlap_unb_ctx = aneuploidy_id %in% overlap_aneuploidy_unb_ctx$aneuploidy_id )



# Export ----
if(length(Sys.glob(chromosome_cn_calls_path))>0 & override_if_exists==F) {
  print("Output already exists, did not override")
} else {
  write.table(cohort_fga_per_chrom,cohort_fga_per_chrom_path,sep="\t",col.names = T,row.names = F)
  write.table(cohort_fga_per_chrom_per_call,cohort_fga_per_chrom_per_call_path,sep="\t",col.names = T,row.names = F)
  write.table(cohort_ploidy,cohort_ploidy_path,sep="\t",col.names = T,row.names = F)
  write.table(complex_cohort_fga_per_chrom,cohort_complex_fga_per_chrom_path,sep="\t",col.names = T,row.names = F)
  
  write.table(chromosome_cn_calls,chromosome_cn_calls_path,sep="\t",row.names = F,col.names = T)
  write.table(overlap_aneuploidy_unb_ctx,overlap_aneuploidy_unb_ctx_path,sep="\t",row.names = F,col.names = T)
}

