#SV hotspots 
#2023-07. 
# Last update 2023-09-13 refactor

#use coverage and peak calling to get recurrently altered regions
# sv bp with window 1mb, also tried 5mb but that worked less well
#then  map back the SVs to the highest peak
#use all svs for making peaks and then when mapping back, distinguish:
#- complex or not
#- cancer type


suppressPackageStartupMessages({
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(dplyr, quietly=TRUE)
  library(GenomicRanges)
  library(openssl)
})



root_dir="~/PycharmProjects/"
source(paste0(root_dir,"structuralvariation/sv_functional_analysis/default.conf"))

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.cna.R"))

source(paste0(wdir,"solids.conf"))

map_template_vars=c('${resources_dir}'=resources_dir,
                    '${merged_svs_dir}'=merged_svs_dir,'${cohort_wdir}'=cohort_wdir,
                    '${cohort_identifier}'=cohort_identifier)



output_dir=paste0(cohort_results_dir,"recurrent_regions/")


#settings  ----
region_size_label="1mb"
window_around_sv=1e6 #1Mb window around SV bp
run_date=20230701

#input
merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.",run_date,".tsv")
merged_svs_path=merged_svs_classes_path
flag_load_from_file=T
#identifiying recurrent region does not require complex sv annotation but currently counting simple/complex separately inside the peaks

#output
override_if_exists=F
#not needed to export
#sv_bp_flank_1mb_coverage_df_path = paste0(output_dir,"sv_bp_",region_size_label,"_peaks.",run_date,".tsv")
#sv_bp_flank_1mb_coverage_complex_df_path = paste0(output_dir,"sv_bp_",region_size_label,"_peaks.complex.",run_date,".tsv") #disabled currently

peak_counts_overview_pancancer_path = paste0(output_dir,"recurrent_regions_svs.counts_overview.pancancer.",region_size_label,".",run_date,".tsv")
peak_counts_overview_cancertype_path = paste0(output_dir,"recurrent_regions_svs.counts_overview.cancertype.",region_size_label,".",run_date,".tsv")
map_svs_highest_peaks_path = paste0(output_dir,"recurrent_regions_svs.map_svs_highest_peaks.",region_size_label,".",run_date,".tsv")


# Cohort ----


cohort = read.table(patient_table_path,sep = "\t", header=T)

patient_id_to_labels = cohort[,c("patient_id","patient_label")]
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")


# Resources ----


## chromosomes -----

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

chromosomes_df = get_chromosomes(chromosome_bands_df)


chrom_order=unique(chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")])



## Genes -----


gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)
gene_properties_df = get_gene_properties_df(gtf_path)

gene_properties_df$gene_coord = gene_properties_df$to_coordinate
gene_properties_df$gene_chrom = gene_properties_df$seqnames
gene_cols=c("gene_id","ensembl_id","gene_name","gene_type","gene_coord","gene_chrom")
gene_cols = gene_cols[gene_cols %in% names(gene_properties_df)]

genes_of_interest_collection = get_genes_of_interest_collection()

gene_properties_df = gene_properties_df %>% 
  mutate(flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
         flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name) 

gene_properties_df_flags = names(gene_properties_df)[grepl("flag",names(gene_properties_df))]
for(flag in gene_properties_df_flags) {
  gene_properties_df[is.na(gene_properties_df[,flag]),flag]=F
}

gene_df_cols = c(gene_cols,gene_properties_df_flags)
gene_df_cols[!gene_df_cols %in% names(gene_properties_df)]


gene_properties=GRanges(gene_properties_df)
names(gene_properties) = gene_properties$gene_id

genes_1mb = GRanges(gene_properties_df)
names(genes_1mb) =genes_1mb$gene_id
genes_1mb =  resize(genes_1mb, width = width(genes_1mb)+1e6, fix = "center")


# Functions -----

get_map_svs_highest_peaks = function(input_svs,target_peaks,input_svs_df=NULL,target_peaks_df=NULL,svs_id_col="patient_sv_merged",peak_id_col="peak_id",target_peak_cols = c("patient_cnt","peak_width"),target_svs_cols = c("patient_label"),cancer_type_split=T) {
  
  target_peak_cols = c(peak_id_col,target_peak_cols) %>% unique()
  target_svs_cols = c(svs_id_col,target_svs_cols) %>% unique()

  if(cancer_type_split) {
    target_peak_cols = c(target_peak_cols,"cancer_type") %>% unique()
    target_svs_cols = c(target_svs_cols,"cancer_type") %>% unique()
  }
  
  if(is.null(input_svs_df)){
    input_svs_df = input_svs %>% as.data.frame() 
  } 
  if(is.null(target_peaks_df)) {
  target_peaks_df = target_peaks %>% as.data.frame()
  }
  
  map_svs_peaks = get_reciprocal_overlap_pairs(input_svs,target_peaks,reciprocal_overlap = 0,svtype_matching = F)
  map_svs_peaks = map_svs_peaks %>% dplyr::rename(!!sym(svs_id_col):=set1,!!sym(peak_id_col):=set2) %>%
    left_join(target_peaks_df[,c(target_peak_cols)],by=peak_id_col) %>% 
    left_join(input_svs_df[,c(target_svs_cols)],by=svs_id_col)
      
  if(cancer_type_split) {
    #match on cancer type depending on the analysis
    map_svs_peaks = map_svs_peaks %>% filter(cancer_type.x==cancer_type.y)
  }
  
  map_svs_highest_peaks = map_svs_peaks %>% group_by(across(all_of(svs_id_col))) %>% slice_max(patient_cnt)
    
    #[map_svs_peaks$patient_cnt == ave(map_svs_peaks$patient_cnt, map_svs_peaks[,svs_id_col],FUN=max),]
  
  duplicates = map_svs_highest_peaks[duplicated(map_svs_highest_peaks$patient_sv_merged),]$patient_sv_merged
  
  resolved_duplicates = map_svs_highest_peaks %>% filter(patient_sv_merged %in% duplicates) %>% group_by(across(all_of(svs_id_col))) %>% slice_max(overlap_set1_set2)
  
  map_svs_highest_peaks = rbind(map_svs_highest_peaks %>% filter(!patient_sv_merged %in% duplicates),resolved_duplicates)
  return(map_svs_highest_peaks)
}

# Read data -----

## Load SV data -----

if(length(Sys.glob(merged_svs_path)) ==1 & flag_load_from_file==T) {
  merged_svs = read.table(merged_svs_path,sep="\t", header=T)
  merged_svs_gr = GRanges(merged_svs$sv_merged_coordinate)
  mcols(merged_svs_gr) = merged_svs %>% select(-start,-end,-seqnames) 
  names(merged_svs_gr) = merged_svs_gr$patient_sv_merged
  
  merged_sv_bp_gr = get_svs_as_bp_gr(merged_svs_gr)

  } else {
  unfiltered_svs_df = load_cohort_svs(cohort,c('${merged_svs_dir}'=merged_svs_dir),svs_path_template=svs_union_anno_somatic_multitool_path_template,flag_add_patient_label = T)
  unfiltered_svs_df = unfiltered_svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
  unfiltered_svs_df = unfiltered_svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
  unfiltered_svs_df = annotate_sv_multitool_support(unfiltered_svs_df,sv_tumor_af_threshold = 0.1)
  unfiltered_svs_df = unfiltered_svs_df %>% unique()
  
  merged_svs = make_merged_svs(unfiltered_svs_df)
  merged_svs = merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  
  
  merged_svs = merged_svs %>% filter(flag_in_filtered_svs) %>% filter(patient_label %in% cohort$patient_label) %>% unique()
  #sets globals therefore no output
  resolve_orphan_svs(merged_svs,unfiltered_svs_df)
  
  #check consistency merged svs and svs df 
  merged_svs %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% nrow() == 0
  unfiltered_svs_df %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% filter(patient_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() == 0
  
  
  merged_svs_classes = read.table(merged_svs_classes_path,sep="\t", header=T)
  merged_svs=merged_svs %>% left_join(merged_svs_classes %>% select(patient_sv_merged,contains("complex")))
  
  
  merged_svs_gr = GRanges(merged_svs$sv_merged_coordinate)
  mcols(merged_svs_gr) = merged_svs 
  names(merged_svs_gr) = merged_svs_gr$patient_sv_merged
  
  merged_svs_coord = GRanges(merged_svs$sv_merged_coordinate) %>% as.data.frame()
  merged_svs_coord$patient_sv_merged = merged_svs$patient_sv_merged
  
  merged_svs = merged_svs %>% left_join(merged_svs_coord[,c("start","end","seqnames","patient_sv_merged")])
            
  merged_svs = merged_svs %>% rowwise() %>% dplyr::mutate(chrom_pair=ifelse(svtype=="CTX",paste0(sort(c(chrom,partner_chrom)),collapse = "-"),NA)) %>% as.data.frame()
  
  merged_sv_bp_gr = get_svs_as_bp_gr(merged_svs_gr)
  
}




#Get recurrent regions by SVs ----

sv_bp_flank_1mb = resize_gr_distance(merged_sv_bp_gr,window_around_sv)
reduced_sv_bp_flank_1mb  = unlist(GenomicRanges::reduce(split(sv_bp_flank_1mb, sv_bp_flank_1mb$patient_id),ignore.strand=T))
reduced_sv_bp_flank_1mb$patient_id=names(reduced_sv_bp_flank_1mb)
sv_bp_flank_1mb_coverage_df = call_peaks(reduced_sv_bp_flank_1mb,max_peaks_only = F)
sv_bp_flank_1mb_coverage_df$peak_type=paste0("sv_bp_",region_size_label)


sv_peaks_overlaps = get_sv_peaks_overlaps(sv_bp_flank_1mb,merged_svs,sv_bp_flank_1mb_coverage_df,svs_id_col="patient_sv_merged",peak_id_col="peak_id")
patients_per_peak = get_patients_per_peak(sv_peaks_overlaps)
sv_bp_flank_1mb_coverage_df = sv_bp_flank_1mb_coverage_df %>% left_join(patients_per_peak)





if(FALSE){

## complex svs only [not used]

sv_bp_flank_1mb_complex = sv_bp_flank_1mb[filter(merged_svs,flag_is_complex)$patient_sv_merged]
reduced_sv_bp_flank_1mb_complex  = unlist(GenomicRanges::reduce(split(sv_bp_flank_1mb_complex, sv_bp_flank_1mb_complex$patient_id),ignore.strand=T))
reduced_sv_bp_flank_1mb_complex$patient_id=names(reduced_sv_bp_flank_1mb_complex)

sv_bp_flank_1mb_coverage_complex_df = call_peaks(reduced_sv_bp_flank_1mb_complex,max_peaks_only = F)

sv_bp_flank_1mb_coverage_complex_df$peak_type=paste0("sv_bp_",region_size_label,"_complex")


sv_peaks_overlaps = get_sv_peaks_overlaps(sv_bp_flank_1mb_complex,merged_svs,sv_bp_flank_1mb_coverage_complex_df,svs_id_col="patient_sv_merged",peak_id_col="peak_id")
patients_per_peak = get_patients_per_peak(sv_peaks_overlaps)
sv_bp_flank_1mb_coverage_complex_df = sv_bp_flank_1mb_coverage_complex_df %>% left_join(patients_per_peak)
  
}


## Reduce SVs separately per cancer type [not used] 
if(FALSE) {
cancer_type_sv_bp_1mb = get_cancer_type_peaks(sv_bp_flank_1mb,peak_type_value=paste0("sv_bp_",region_size_label))
#cancer_type_sv_bp_1mb_complex = get_cancer_type_peaks(sv_bp_flank_1mb_complex,peak_type_value=paste0("sv_bp_",region_size_label,"_complex")
}


# Get recurrent regions - map svs to highest peak ----
#map to highest with ave and then summarize on peak level for better estimate of svs cnt
#hope to get single peaks that way instead of staircase


peak_id_col="peak_id"
sv_bp_1mb_peaks_gr= GRanges(sv_bp_flank_1mb_coverage_df %>% select(-patient_md5))
names(sv_bp_1mb_peaks_gr) = sv_bp_1mb_peaks_gr$peak_id

#to make sure not duplicated merged_sv_bp_gr[merged_svs$patient_sv_merged]
map_svs_highest_peaks = get_map_svs_highest_peaks(input_svs = merged_sv_bp_gr[merged_svs$patient_sv_merged],
                                                  target_peaks = sv_bp_1mb_peaks_gr,
                                                  input_svs_df = merged_svs,
                                                  target_svs_cols = c("patient_label","cancer_type","complex_sv_id","complex_sv_class_super","flag_is_complex"),
                                                  cancer_type_split = F)


gene_df_cols = c("gene_id","gene_name","gene_coord","gene_type","seqnames",gene_properties_df_flags)

map_peaks_genes = get_reciprocal_overlap_pairs(sv_bp_1mb_peaks_gr,gene_properties,reciprocal_overlap = 0,svtype_matching = F)  %>% 
  dplyr::rename(peak_id=set1,gene_id=set2) %>% 
  left_join(gene_properties_df[,gene_df_cols] %>% 
              dplyr::rename(gene_chrom=seqnames))


map_peaks_genes_1mb = get_reciprocal_overlap_pairs(sv_bp_1mb_peaks_gr,genes_1mb,reciprocal_overlap = 0,svtype_matching = F)  %>%
  dplyr::rename(peak_id=set1,gene_id=set2) %>% 
  left_join(gene_properties_df[,gene_df_cols] %>% 
              dplyr::rename(gene_chrom=seqnames))



peaks_all_svs = map_svs_highest_peaks %>% get_patients_per_peak()
peaks_complex_svs = map_svs_highest_peaks %>% get_patients_per_peak(group_cols = c(peak_id_col,"flag_is_complex"))
peaks_svs_cancer_type = map_svs_highest_peaks %>% get_patients_per_peak(group_cols = c(peak_id_col,"cancer_type"),no_cancer_type=T)
peaks_complex_svs_cancer_type = map_svs_highest_peaks %>% get_patients_per_peak(group_cols = c(peak_id_col,"flag_is_complex","cancer_type"),no_cancer_type=T)
peaks_complex_class_svs_cancer_type = map_svs_highest_peaks %>% get_patients_per_peak(group_cols = c(peak_id_col,"flag_is_complex","complex_sv_class_super","cancer_type"),no_cancer_type=T)

# Make overview of hotspots   ----
#overview starting from all peaks list (pan cancer or split per cancer type) and then join additional cols such that only 1 row per peak 

peak_coord = sv_bp_flank_1mb_coverage_df %>% select(peak_id,seqnames,start,end)

peak_counts_overview_pancancer = peaks_all_svs %>% 
  dplyr::rename(all = patient_cnt, all_cancer_type=cancer_type_cnt) %>% 
  left_join(peaks_complex_svs %>% pivot_wider(id_cols="peak_id",names_from = "flag_is_complex",values_from = "patient_cnt",values_fill=0) %>% dplyr::rename(simple=`FALSE`,complex=`TRUE`)) %>%
  left_join(peaks_complex_svs %>% pivot_wider(id_cols="peak_id",names_from = "flag_is_complex",values_from = "cancer_type_cnt",values_fill=0) %>% dplyr::rename(simple_cancer_type=`FALSE`,complex_cancer_type=`TRUE`)) %>%
  left_join(peaks_complex_svs %>% filter(flag_is_complex) %>% dplyr::rename(complex_patient_lst=patient_lst) %>% select(peak_id,complex_patient_lst)) %>%
  left_join(peak_coord) %>%
  left_join(map_peaks_genes %>% filter(flag_pmc_panel) %>% group_by(peak_id) %>% summarize(genes_lst_pmc=lst_str(gene_name))) %>%
  left_join(map_peaks_genes_1mb %>% filter(flag_pmc_panel) %>% group_by(peak_id) %>% summarize(genes_1mb_lst_pmc=lst_str(gene_name))) %>%
  left_join(map_peaks_genes %>% filter(flag_gene_of_interest) %>% group_by(peak_id) %>% summarize(genes_lst=lst_str(gene_name))) %>%
  left_join(map_peaks_genes_1mb %>% filter(flag_gene_of_interest) %>% group_by(peak_id) %>% summarize(genes_1mb_lst=lst_str(gene_name))) 

top_recurrent_peaks_pancancer = peak_counts_overview_pancancer %>% filter(complex>2) 


peak_counts_overview_cancertype = peaks_svs_cancer_type %>% 
  dplyr::rename(all = patient_cnt) %>% 
  left_join(peaks_complex_svs_cancer_type %>% pivot_wider(id_cols=c("peak_id","cancer_type"),names_from = "flag_is_complex",values_from = "patient_cnt",values_fill=0) %>% dplyr::rename(simple=`FALSE`,complex=`TRUE`)) %>%
  left_join(peaks_complex_svs_cancer_type %>% filter(flag_is_complex) %>% dplyr::rename(complex_patient_lst=patient_lst) %>% select(peak_id,complex_patient_lst,cancer_type)) %>%
  left_join(peak_counts_overview_pancancer %>% select(peak_id,all,complex,simple) %>% dplyr::rename_with(.cols=-peak_id,function(x){paste0("pancancer_",x)})) %>%
  left_join(peak_coord) %>% 
  left_join(map_peaks_genes %>% filter(flag_pmc_panel) %>% group_by(peak_id) %>% summarize(genes_lst_pmc=lst_str(gene_name))) %>%
  left_join(map_peaks_genes_1mb %>% filter(flag_pmc_panel) %>% group_by(peak_id) %>% summarize(genes_1mb_lst_pmc=lst_str(gene_name))) %>%
  left_join(map_peaks_genes %>% filter(flag_gene_of_interest) %>% group_by(peak_id) %>% summarize(genes_lst=lst_str(gene_name))) %>%
  left_join(map_peaks_genes_1mb %>% filter(flag_gene_of_interest) %>% group_by(peak_id) %>% summarize(genes_1mb_lst=lst_str(gene_name))) 

top_recurrent_peaks_cancertype = peak_counts_overview_cancertype %>% filter(complex>2) 



## export ----

if(length(Sys.glob(peak_counts_overview_pancancer_path))>0 & override_if_exists==F) {
  print("Output already exists, did not override")
} else {
#no need to export these
#write.table(sv_bp_flank_1mb_coverage_df,sv_bp_flank_1mb_coverage_df_path,quote=F,row.names = F,col.names = T,sep="\t")
#write.table(sv_bp_flank_1mb_coverage_complex_df,sv_bp_flank_1mb_coverage_complex_df_path,quote=F,row.names = F,col.names = T,sep="\t")

write.table(peak_counts_overview_pancancer,peak_counts_overview_pancancer_path,sep = "\t",col.names = T,row.names = F,quote=F)
write.table(peak_counts_overview_cancertype,peak_counts_overview_cancertype_path,sep = "\t",col.names = T,row.names = F,quote=F)
write.table(map_svs_highest_peaks,map_svs_highest_peaks_path,sep = "\t",col.names = T,row.names = F,quote=F)
}

