## Get amplicons 
## Last update: 2023-06-06 
# make amplicons for cohort
# based on copy ratio peaks relative to copy ratio mean of chromosome
# annotate with genes
# export merged seg and footprints + genes
# => pediatric cancer genes, gene of interest collection, PMC panel (maxima list)


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
#defaults?

#in solids conf
#cr_amplification_threshold=1.9
#autosomes= c(paste("chr",1:22,sep=""),"chrX")

amplicons_footprint_size_min = 50e3 # min size of amplicon to be included in analysis

## Cohort paths ----


#complex_svs_results_dir=paste0(cohort_results_dir,"complex_sv_clusters/") #todo replace
       
## Outputs

rundate=20230701
amplicon_merged_segments_path = paste0(cohort_results_dir,"amplicon_merged_segments.",rundate,".tsv")
amplicon_footprints_path = paste0(cohort_results_dir,"amplicon_footprints.",rundate,".tsv")

# Cohort ----

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

automated_baseline_shifts = read.table(automated_baseline_shifts_path,sep="\t",header=T)
manual_baseline_shifts = read.table(manual_baseline_shifts_path,sep="\t",header=T)
baseline_shifts= get_baseline_correction(cohort,automated_baseline_shifts,manual_baseline_shifts)
baseline_shifts$flag_is_diploid_no_relative_correction=NA

# RESOURCES----

## Chromosomes ----
map_template_vars=c('${resources_dir}'=resources_dir)
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)
chr_arms_df = get_chr_arms(chromosome_bands_df)
chromosomes_df = get_chromosomes(chr_arms_df)
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames
#chr_arm_order = chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")]# %>% flatten_chr()
#chrom_order = chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")] %>% unique()# %>% flatten_chr()

### Load genes ----

gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)
gene_properties_df = get_gene_properties_df(gtf_path)

gene_properties_df$gene_coord = gene_properties_df$to_coordinate
gene_properties_df$gene_chrom = gene_properties_df$seqnames
gene_cols=c("gene_id","ensembl_id","gene_name","gene_type","gene_coord","gene_chrom")
gene_cols = gene_cols[gene_cols %in% names(gene_properties_df)]

gene_properties=GRanges(gene_properties_df)
names(gene_properties) = gene_properties$gene_id

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



# Prepare data ----
## Load CN data ----
segments_df = load_cohort_cn_segments(cohort,map_template_vars=map_template_vars,baseline_shifts=baseline_shifts,flag_add_patient_label=T)
segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id

chrom_centered_segments_df = get_chrom_centered_cn_segments(segments,segments_df,chromosomes,chromosomes_df,cohort) 
chrom_centered_segments = GRanges(chrom_centered_segments_df)
names(chrom_centered_segments) = chrom_centered_segments$patient_cna_id

## Make amplicons ----


amplicon_segments = chrom_centered_segments[chrom_centered_segments$cr_l2fc_50>cr_amplification_threshold & chrom_centered_segments$cr_l2fc_50_uncorrected>cr_amplification_threshold & chrom_centered_segments@seqnames %in% autosomes]
# also uncorrected >2 otherwise in deep deletion "amp"
#chrom_centered_segments_df %>% filter(patient_label=="M763AAB"  & seqnames=="chr4" & cr_l2fc_50>2) %>% select(patient_cna_id,cr_l2fc_50,chrom_mean_cr)
#not selected due to additional criterium
#chrom_centered_segments_df %>% filter(cr_l2fc_50_uncorrected > 2 & cr_l2fc_50<2) %>% select(patient_label,seqnames,cr_l2fc_50,chrom_mean_cr,tumor_id,)
  

amplicon_merged_seg = get_cna_merged_gainloss_per_patient(amplicon_segments,chrom_centered_segments_df,
                                                           cnas_df_cols = c("patient_cna_id","cr_l2fc_50","maf_50","call","cna_coordinate","chrom_mean_cr"),
                                                           group_cols = c("chrom_mean_cr"))
amplicon_merged_seg$cr_l2fc_50 = amplicon_merged_seg$cr_l2fc_50+amplicon_merged_seg$chrom_mean_cr
amplicon_merged_seg$cr_l2fc_50_max = amplicon_merged_seg$cr_l2fc_50_max+amplicon_merged_seg$chrom_mean_cr

amplicon_merged_seg_df = amplicon_merged_seg %>% as.data.frame()
amplicon_merged_seg_df$patient_cna_merged_id = rownames(amplicon_merged_seg_df)
amplicon_merged_seg_df$chrom = as.character(amplicon_merged_seg_df$seqnames)
amplicon_merged_seg_df$footprint_id=paste0(amplicon_merged_seg_df$patient_label,"_",amplicon_merged_seg_df$chrom)
amplicon_merged_seg$patient_cna_merged_id=names(amplicon_merged_seg)
mcols(amplicon_merged_seg) = mcols(amplicon_merged_seg) %>% as.data.frame() %>% left_join(amplicon_merged_seg_df %>% select(-seqnames,-start,-end,-width,-strand))

amplicon_merged_seg_df = amplicon_merged_seg_df  %>% left_join(cohort[,c("patient_label","cancer_type")])


amplicon_footprints = amplicon_merged_seg_df %>%
  rowwise() %>%
  group_by(footprint_id,patient_label,seqnames,cancer_type,chrom_mean_cr) %>%
  summarize(max_cr=max(cr_l2fc_50_max),
            median_cr=median(cr_l2fc_50),
            peaks=cnt_str(patient_cna_merged_id),
            seqnames=unique(seqnames),min_start=min(start),max_end=max(end),
            footprint_size=max_end-min_start,
            footprint_coordinate=paste0(seqnames,":",min_start,"-",max_end)) 

amplicon_footprints$footprint_width=amplicon_footprints$footprint_size

amplicon_footprints = amplicon_footprints %>% filter(footprint_size>amplicons_footprint_size_min)

#remove very small => I checked and they do not overlap cancer genes
amplicon_merged_seg_df = amplicon_merged_seg_df %>% filter(footprint_id %in% amplicon_footprints$footprint_id)



## Amplicon overlap with Genes ---
#amplicon_merged_seg also contains small gaps so use proper segments in this region

map_amplicons_genes = get_reciprocal_overlap_pairs(amplicon_merged_seg,gene_properties,reciprocal_overlap = 0,svtype_matching = F)
map_amplicons_genes = map_amplicons_genes %>% dplyr::rename(patient_cna_merged_id=set1,gene_id=set2) 

genes_cn_values_df = get_chrom_cn_fractions(segments,segments_df,
                                            gene_properties[gene_properties$gene_id %in% unique(map_amplicons_genes$gene_id)],
                                            gene_properties_df,cohort,ranges_id_col="gene_id",ranges_width_col="gene_width",
                                            relative_cr = F,return_wide = T)

genes_cn_values_df = genes_cn_values_df %>% select(patient_label,gene_id,call,cr_l2fc_50) %>% unique()
                                                                                          
map_amplicons_genes = map_amplicons_genes  %>% 
  left_join(amplicon_merged_seg_df %>% select(patient_label,patient_cna_merged_id,footprint_id,patient_label,cancer_type,seqnames),by="patient_cna_merged_id") %>%
  left_join(genes_cn_values_df) %>%
  left_join(gene_properties_df[,gene_df_cols],by="gene_id") 

amplicons_genes_summary = map_amplicons_genes %>% 
  filter(flag_gene_of_interest) %>%
  group_by(footprint_id,patient_label,cancer_type,seqnames) %>% 
  summarize(genes_lst= toString(sort(unique(paste0(gene_name," (",call," cr:",round(cr_l2fc_50,2),", frac:",round(overlap_set2_set1,2),")")))))


# Export ----
amplicon_footprints = amplicon_footprints %>% left_join(amplicons_genes_summary) 

write.table(amplicon_footprints,amplicon_footprints_path,sep = "\t",col.names = T,row.names = F)
write.table(amplicon_merged_seg_df,amplicon_merged_segments_path,sep = "\t",col.names = T,row.names = F)



