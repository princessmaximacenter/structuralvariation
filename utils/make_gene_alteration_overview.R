

wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))

source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(script_dir,"functions.cna.R"))
source(paste0(script_dir,"functions.expression.R"))

#source(paste0(utils_script_dir,"functions.cohort_analyses.R"))

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


## Fill path templates
map_template_vars=c('${resources_dir}'=resources_dir,'${merged_svs_dir}'=merged_svs_dir,'${utils_output_dir}'=utils_output_dir,
                    '${cohort_wdir}'=cohort_wdir,'${cohort_identifier}'=cohort_identifier)

source(paste0(wdir,"solids.conf"))

sv_tumor_af_threshold=0.1
snv_tumor_af_threshold=0.1
autosomes= c(paste("chr",1:22,sep=""),"chrX")
vst_zscore_threshold=1.98
snv_promoter_hotspot_lst = c("chr5:1295113_G/A","chr5:1295135_G/A")


l2fc_wilcox_threshold = log2(1.1)
fdr_wilcox_threshold = 0.2

#gene_centric_results_dir= "~/results/gene_centric/"
gene_centric_results_dir = "~/data/sv/"

cohort_results_dir
# todo make results here:




merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.20230701.tsv")


#input
amplicon_footprints_path = paste0(cohort_results_dir,"amplicon_footprints.20230701.tsv")
amplicon_merged_segments_path = paste0(cohort_results_dir,"amplicon_merged_segments.20230701.tsv")


run_date=20230701
cohort_snv_path = paste0(cohort_results_dir,"cohort_snv.somatic.tsv")
cohort_snv_burden_somatic_path =paste0(cohort_results_dir,"cohort_snv_burden_somatic.",run_date,".tsv")

# output
override_if_exists=F

patient_gene_overview_detailed_path = paste0(cohort_results_dir,"gene_level_overview.somatic.detailed.tsv")
patient_gene_overview_affected_path = paste0(cohort_results_dir,"gene_level_overview.somatic.affected.tsv")
patient_gene_overview_snv_sv_path = paste0(cohort_results_dir,"gene_level_overview.somatic.snv_sv_bp.tsv")
patient_gene_overview_no_impact_filter_path = paste0(cohort_results_dir,"gene_level_overview.somatic.snv_sv_bp.no_impact_filter.tsv")

patient_gene_overview_detailed_cancer_path = paste0(cohort_results_dir,"gene_level_overview.somatic.detailed.cancer.tsv")



#settings 
cr_hom_loss_threshold = -1.5 #-1 is 1 copy left so less than 1 copy => 0.7, use abs cr so not ploidy corr
cr_amplification_threshold=1.9 #same as amplicon analysis
cr_minimum_cn_change=0.2

# Resources ----
## Functions ----

get_gene_cna_from_segments = function(segments,segments_df,gene_properties,gene_properties_df) {
  gene_cna_overlaps = get_reciprocal_overlap_pairs(segments,gene_properties,reciprocal_overlap = 0,svtype_matching = F)
  gene_cna_overlaps = gene_cna_overlaps %>% dplyr::rename(gene_id=set2, patient_cna_id=set1)
  
  cnas_df_cols = c("patient_label","patient_cna_id","cr_l2fc_50","maf_50","call")
  gene_cols=c("gene_id","ensembl_id","gene_name","gene_type")
  
  gene_cna_overlaps = gene_cna_overlaps %>% 
    left_join(segments_df[,cnas_df_cols]) %>% 
    left_join(gene_properties_df[,gene_cols])
  
  # summarize over the feature
  gene_cna_overlaps = gene_cna_overlaps %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4)
  
  group_cols=c("patient_label","gene_name","gene_id","gene_type","ensembl_id")
  gene_cna = get_overlaps_modeled_seg_summary_short(gene_cna_overlaps,group_cols)
  gene_cna = gene_cna %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4)
  return(gene_cna)
}

## Cancer genes ----

cancer_genes = get_cancer_genes(resources_dir)
cancer_related=get_cancer_genes(resources_dir) %>% dplyr::rename(gene_name=gene_id)

genes_of_interest = cancer_related[,c("gene_name")] %>% as.data.frame()
names(genes_of_interest) = c("gene_name")

gene_onco_tsg_consensus = get_gene_onco_tsg_consensus(cancer_genes)




## Chromosome bands and arms ----

chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df=as.data.frame(chr_arms)

chromosome_bands_df = chromosome_bands_df %>% mutate(chr_arm = ifelse(grepl("p",cytoband),paste0(seqnames,"p"),paste0(seqnames,"q"))) 

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband

chromosomes_df = get_chromosomes(chr_arms_df)
chromosomes = GRanges(chromosomes_df)
names(chromosomes)=chromosomes_df$seqnames


## Genes ----
gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)

gene_properties_df = get_gene_properties_df(gtf_path)

gene_properties=GRanges(gene_properties_df)
names(gene_properties) = gene_properties$gene_id

genes_to_cytobands = get_genes_to_cytobands(gene_properties,chromosome_bands)


gene_properties_df = gene_properties_df %>% left_join(genes_to_cytobands)
gene_properties_df$gene_coord = gene_properties_df$to_coordinate
gene_properties_df$gene_name_display =  paste0(gene_properties_df$gene_name,"_",gene_properties_df$cytoband,"_",gene_properties_df$ensembl_id)
gene_cols=c("gene_id","ensembl_id","gene_name","gene_type","gene_coord","seqnames","cytoband")
gene_cols = gene_cols[gene_cols %in% names(gene_properties_df)]

genes_of_interest_collection = get_genes_of_interest_collection()
 
gene_properties_df = gene_properties_df %>% 
  mutate(flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
         flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name) 

gene_properties_df = gene_properties_df %>%  left_join(cancer_genes[,c("gene_id","oncogene","tsg")],by=c("ensembl_id"="gene_id")) %>% dplyr::rename(flag_oncogene=oncogene,flag_tsg=tsg)

#check gene_properties_df[duplicated(gene_properties_df$gene_id),]

gene_properties_df_flags = names(gene_properties_df)[grepl("flag",names(gene_properties_df))]

for(flag in gene_properties_df_flags) {
  gene_properties_df[is.na(gene_properties_df[,flag]),flag]=F
}

gene_df_cols = c(gene_cols,gene_properties_df_flags)
gene_df_cols[!gene_df_cols %in% names(gene_properties_df)]


genes_1mb = GRanges(gene_properties_df)
names(genes_1mb) =genes_1mb$gene_id
genes_1mb =  resize(genes_1mb, width = width(genes_1mb)+1e6, fix = "center")


# Cohort ----
##Load cohort and annotate 

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)

patient_id_to_labels = cohort[,c("patient_id","patient_label")]

patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

cohort_anno_display_cols=patient_tumor_id_cols

automated_baseline_shifts = read.table(automated_baseline_shifts_path,sep="\t",header=T)
manual_baseline_shifts = read.table(manual_baseline_shifts_path,sep="\t",header=T)
baseline_shifts= get_baseline_correction(cohort,automated_baseline_shifts,manual_baseline_shifts)

# Load Cohort SVs ----

merged_svs = read.table(merged_svs_classes_path,sep="\t",header=T)

merged_svs_gr = GRanges(merged_svs$sv_merged_coordinate)
mcols(merged_svs_gr) = merged_svs %>% select(-start,-end,-seqnames) 
names(merged_svs_gr) = merged_svs_gr$patient_sv_merged

merged_sv_bp_gr = get_svs_as_bp_gr(merged_svs_gr)


## overlap svs with genes ----
#sanity check all with sv also in our sv merged anno

overlaps_genes_svs = get_reciprocal_overlap_pairs(gene_properties,merged_sv_bp_gr,ignore_strand = T,svtype_matching = F,reciprocal_overlap=0) %>%
  dplyr::rename(gene_id=set1, patient_sv_merged = set2 ) %>% 
  left_join(gene_properties_df[,gene_df_cols]) %>%
  left_join(merged_svs %>% dplyr::select(patient_sv_merged,"svtype","tumor_af_mean","svlen_mean_str",contains("complex"),patient_label))

genes_sv_bp = overlaps_genes_svs %>%
  group_by(patient_label,gene_name,ensembl_id) %>%
  summarize(sv_bp = toString(paste0(svtype," ",round(tumor_af_mean,2)," taf ",svlen_mean_str)),
            sv_bp_cnt = n(),
            sv_bp_max_taf=max(tumor_af_mean))


#overlaps_genes_svs_1mb = get_reciprocal_overlap_pairs_start_end(merged_svs_gr,genes_1mb,ignore_strand = T,svtype_matching = F,reciprocal_overlap=0) %>%
overlaps_genes_svs_1mb = get_reciprocal_overlap_pairs(genes_1mb,merged_sv_bp_gr,ignore_strand = T,svtype_matching = F,reciprocal_overlap=0) %>%
  dplyr::rename(gene_id=set1, patient_sv_merged = set2 ) %>% 
  left_join(gene_properties_df[,gene_df_cols]) %>%
  left_join(merged_svs %>% dplyr::select(patient_sv_merged,svtype,tumor_af_mean,svlen_mean_str,contains("complex"),patient_label,sv_merged_coordinate))

overlaps_genes_svs_1mb$distance = distance(GRanges(overlaps_genes_svs_1mb$gene_coord),GRanges(overlaps_genes_svs_1mb$sv_merged_coordinate))

#2023-08-22 update subtract patient_sv_merged
genes_sv_1mb = overlaps_genes_svs_1mb %>% 
  anti_join(overlaps_genes_svs[,c("patient_sv_merged","gene_id")]) %>%
  group_by(patient_label,gene_name,ensembl_id) %>%
  summarize(sv_1mb = toString(paste0(svtype," ",round(tumor_af_mean,2)," taf dist.",round(distance/1e3,2)," kbp")),
          sv_1mb_cnt = n(),
          sv_1mb_max_taf=max(tumor_af_mean))

#2023-08-22 add complex svs
complex_genes_svs = overlaps_genes_svs %>% 
  filter(flag_is_complex) %>%
  group_by(patient_label,gene_name,ensembl_id) %>%
  summarize(complex_sv_bp = toString(paste0(svtype," ",round(tumor_af_mean,2))),
            complex_sv_bp_cnt = n(),
            complex_sv_bp_max_taf = max(tumor_af_mean))


complex_genes_svs_1mb = overlaps_genes_svs_1mb %>% 
  anti_join(overlaps_genes_svs[,c("patient_sv_merged","gene_id")]) %>%
  filter(flag_is_complex) %>%
  group_by(patient_label,gene_name,ensembl_id) %>%
  summarize(complex_sv_1mb =toString(paste0(svtype," ",round(tumor_af_mean,2)," taf dist.",round(distance/1e3,2)," kbp")),
            complex_sv_1mb_cnt = n(),
            complex_sv_1mb_max_taf = max(tumor_af_mean))


#can only be 1 because of 1mb?
complex_genes_class = overlaps_genes_svs_1mb %>% 
  filter(flag_is_complex) %>%
  group_by(patient_label,gene_name,ensembl_id) %>%
  summarize(complex_sv_class=lst_str(complex_sv_class),
            complex_sv_id=lst_str(complex_sv_id),
            complex_sv_cnt=cnt_str(complex_sv_id))
complex_genes_class %>% filter(complex_sv_cnt>1) %>% nrow() ==0



# somatic SNVs ----

#2022-10-22 refactored but slow
#if(FALSE) {
cohort_snv_somatic = read.table(cohort_snv_path,sep="\t",header=T)
cohort_snv_somatic = cohort_snv_somatic %>% left_join(patient_id_to_labels)
  cohort_snv_somatic = cohort_snv_somatic %>% filter(patient_label %in% cohort$patient_label)
  cohort_snv_somatic$patient_snv_id=paste0(cohort_snv_somatic$patient_label,"_",cohort_snv_somatic$snv_id)
#}

if(FALSE) {
cohort_snv_somatic = load_cohort_snvs(cohort,map_template_vars = c('${utils_output_dir}'=snvs_dir,'${analysis_type}'="somatic"),
                                 analysis_type="somatic",snv_path_template=snv_output_path_template,flag_add_patient_label=T)

#write.table(cohort_snv_somatic,cohort_snv_path,quote=F,sep="\t",row.names = F,col.names = T)

}

cohort_snv_somatic = cohort_snv_somatic %>% mutate(dbsnp_present = grepl("rs",Existing_variation))
cohort_snv_somatic = cohort_snv_somatic %>% mutate(cosmic_present = grepl("COSV|COSM",Existing_variation)) #legacy identifier COSM added

cohort_snv_somatic = cohort_snv_somatic %>% 
  mutate(flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
         flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name) 

cohort_snv_somatic = cohort_snv_somatic %>% separate(col=PolyPhen,sep="\\(",into=c("polyphen_label","score"),remove=F)
cohort_snv_somatic = cohort_snv_somatic %>% separate(col=SIFT,sep="\\(",into=c("sift_label","score"),remove=F)

cohort_snv_somatic = cohort_snv_somatic %>% mutate(flag_protein_function_affected = ( grepl("damaging",PolyPhen) | grepl("deleterious",SIFT)))
cohort_snv_somatic = cohort_snv_somatic %>% mutate(flag_benign = ( grepl("benign",PolyPhen) | grepl("tolerated",SIFT)))

#always filter by tumor AF, autosomes+chrX, remove blacklisted , remove dbsnp if not in cosmic as well 
#cohort_snv_somatic_unfiltered = cohort_snv_somatic #cohort_snv_somatic_bk in other file

## note 2023-06-29 got a lot more tp53 mutations with new annotation more cosmic ids added that override the dbsnp rs ids 
#consider dropping the filtering on dbsnp  presence

if(cohort_snv_somatic$patient_label %>% unique() %>% length() != nrow(cohort)) {
  print("Prefiltering already missing patients")
}

# some older versions did not include this
# TODO rerun snvs from fusion sq => done
if(cohort_snv_somatic[is.na(cohort_snv_somatic$dac_blacklist),] %>% nrow() > 0) {
  
  blacklist_path="~/Documents/resources/hg38-blacklist.v2.bed"
  
  blacklist = rtracklayer::import.bed(blacklist_path)
  names(blacklist)=paste0("blacklist_region_",1:length(blacklist))
  
  snv_gr = GRanges(cohort_snv_somatic[,c("start","end","seqnames","patient_snv_id")])
  names(snv_gr) = snv_gr$patient_snv_id
  snv_in_blacklist = subsetByOverlaps(snv_gr,blacklist)
  
  cohort_snv_somatic = cohort_snv_somatic %>% mutate(dac_blacklist= patient_snv_id %in% snv_in_blacklist$patient_snv_id)
}


cohort_snv_somatic = cohort_snv_somatic %>% filter(tumor_AF>snv_tumor_af_threshold)
cohort_snv_somatic = cohort_snv_somatic %>% filter(seqnames %in% autosomes) %>% filter(!dac_blacklist)

#for reporting mutations in genes ignore dbsnp filtering
unfiltered_cohort_snv_somatic = cohort_snv_somatic

#remove dbsnp unless also cosmic +> disabled as of 2023-06-29
#cohort_snv_somatic = cohort_snv_somatic %>% filter(!dbsnp_present | cosmic_present)

if(cohort_snv_somatic$patient_label %>% unique() %>% length() != nrow(cohort)) {
  print("Patients with 0 SNV/indels")
  patients_no_snvs = cohort %>% filter(!patient_label %in% cohort_snv_somatic$patient_label)
}
### DIVERGENCE OF GENE AFFECTING AND BURDEN
cohort_snv_somatic_first_filtering= cohort_snv_somatic 

# filtering for effect on gene 
#Filter by impact for if affects gene
cohort_snv_somatic = cohort_snv_somatic_first_filtering %>% filter(IMPACT %in% c("MODERATE","HIGH"))

#remove benign unless also damaging (conflict sift and polyphen) or cosmic
cohort_snv_somatic = cohort_snv_somatic %>% filter(!flag_benign | flag_protein_function_affected | cosmic_present)

## TODO 2023-08-22 after update unnest listed CSQ data get many more rows for single snv_id
# cnt unique snv ids 
# make sure not single snv id ends up  mutating >1 gene of interest.

gene_snv_tbl = cohort_snv_somatic %>% filter(!is.na(gene_name) & gene_name != "" & flag_gene_of_interest) %>% select(gene_name,patient_snv_id) %>% unique()
gene_snv_tbl[duplicated(gene_snv_tbl$patient_snv_id),] %>% nrow() == 0
#is okay

cohort_snv_somatic_filtered_path =paste0(cohort_results_dir,"cohort_snv_somatic.filtered.tsv")
#write.table(cohort_snv_somatic,cohort_snv_somatic_filtered_path,quote=F,sep="\t",row.names = F,col.names = T)

##2023-04-06 add hotspot promoter mutations => TERT
unfiltered_cohort_snv_somatic = unfiltered_cohort_snv_somatic %>%
  dplyr::mutate(polyphen_label=ifelse(snv_id %in% snv_promoter_hotspot_lst,"probably_damaging",polyphen_label),
                sift_label=ifelse(snv_id %in% snv_promoter_hotspot_lst,"deleterious",sift_label))



## SNV burden ----


if(length(Sys.glob(cohort_snv_burden_somatic_path))==1) {
  cohort_snv_burden_somatic = read.table(cohort_snv_burden_somatic_path,sep="\t",header=T)
  #manual 
  cohort_snv_burden_somatic %>% summary()
  cohort_snv_burden_somatic$snv_indel_somatic_nonsyn_cnt %>% summary()
  cohort_snv_burden_somatic$snv_somatic_nonsyn_cnt %>% summary()
  
} else {
  vep104 = rtracklayer::import(paste0(resources_dir,"Homo_sapiens.GRCh38.104.gtf.gz"))
  genes_vep104 = vep104[vep104$type=="gene"]
  seqlevelsStyle(genes_vep104)="UCSC"
  genes_vep104_df = genes_vep104 %>% as.data.frame()
  genes_vep104_df$ensembl_id=genes_vep104_df$gene_id

  # denominator is coding sequence
denominator = 41.072372 ##  changed to 41.072372 because value in pipeline > TODO rerun

## checked what would be removed and decided to not filter burden by DP/AD only for the mutsig calling 
## filter DP tumor and normal >20, remove alt_supporting_read_in_normal
#filter_out_questionmark = cohort_snv_somatic %>% filter(tumor_DP<20|normal_DP<20 | !grepl(", 0)",normal_AD))

burden_cohort_snv_somatic = cohort_snv_somatic_first_filtering %>% filter(seqnames !="chrX")

## choose to remove LOW in its entirety because the splice affecting mutations are not in coding REGION 
burden_cohort_snv_somatic_nonsynonymous = burden_cohort_snv_somatic  %>% 
  filter(IMPACT %in% c("HIGH","MODERATE") & !Consequence %in% c("regulatory_region_ablation"))

burden_cohort_snv_somatic_nonsynonymous = burden_cohort_snv_somatic_nonsynonymous %>% left_join(genes_vep104_df[,c("ensembl_id","gene_biotype")])
burden_cohort_snv_somatic_nonsynonymous = burden_cohort_snv_somatic_nonsynonymous %>% filter( gene_biotype =="protein_coding")

burden_cohort_snv_somatic_nonsynonymous_only_snv =  burden_cohort_snv_somatic_nonsynonymous %>% filter(VARIANT_CLASS=="SNV")

cohort_snv_burden_somatic = burden_cohort_snv_somatic %>% cnt_uq(uq_id = "patient_snv_id", attr_name = "snv_indel_somatic_cnt") 
cohort_snv_burden_somatic = cohort_snv_burden_somatic %>% left_join( burden_cohort_snv_somatic_nonsynonymous %>% cnt_uq(uq_id = "patient_snv_id", attr_name = "snv_indel_somatic_nonsyn_cnt") ) 
cohort_snv_burden_somatic = cohort_snv_burden_somatic %>% left_join( burden_cohort_snv_somatic_nonsynonymous_only_snv %>% cnt_uq(uq_id = "patient_snv_id", attr_name = "snv_somatic_nonsyn_cnt") ) 
cohort_snv_burden_somatic[is.na(cohort_snv_burden_somatic)]=0


cohort_snv_burden_somatic$tmb = cohort_snv_burden_somatic$snv_indel_somatic_nonsyn_cnt/denominator


cohort_snv_burden_somatic = cohort_snv_burden_somatic %>% left_join(patient_id_to_labels)
cohort %>% filter(!patient_label %in% cohort_snv_burden_somatic$patient_label) %>% nrow() == 0

write.table(cohort_snv_burden_somatic,cohort_snv_burden_somatic_path,quote=F,sep="\t",row.names = F,col.names = T)
}



# Load CN segments ----
segments_df = load_cohort_cn_segments(cohort,map_template_vars=map_template_vars,baseline_shifts=baseline_shifts,flag_add_patient_label=T)
segments = GRanges(segments_df)
names(segments) = segments$patient_cna_id

chrom_cn_fractions = get_chrom_cn_fractions(segments,segments_df,
                                            chromosomes,chromosomes_df,cohort,ranges_id_col="seqnames",ranges_width_col="chrom_width",
                                            relative_cr = F,return_wide = T)
chrom_cn_fractions = chrom_cn_fractions %>% dplyr::mutate(chrom_mean_cr=cr_l2fc_50)  
chrom_cn_fractions = chrom_cn_fractions %>% select(patient_label,seqnames,chrom_mean_cr) %>% ungroup()


## amplicons ----
#amplicon_footprints >50 kb 

amplicon_footprints = read.table(amplicon_footprints_path,sep = "\t",header=T)
amplicon_footprints = amplicon_footprints %>% filter(footprint_size>50e3)

amplicon_merged_segments_df = read.table(amplicon_merged_segments_path,sep = "\t",header=T)
amplicon_merged_segments_df = amplicon_merged_segments_df %>% filter(footprint_id %in% amplicon_footprints$footprint_id)

amplicon_merged_seg = GRanges(amplicon_merged_segments_df$cna_coordinate)
#mcols(amplicon_merged_seg) = amplicon_merged_segments_df %>% select(-seqnames,-start,-end)
names(amplicon_merged_seg) = amplicon_merged_segments_df$patient_cna_merged_id
amplicon_merged_segments_df = amplicon_merged_segments_df %>% dplyr::rename_with(.cols=c("cr_l2fc_50","footprint_id","cna_width","cna_coordinate"),.fn=function(x){paste0("amplicon_",x)})

map_amplicons_genes = get_reciprocal_overlap_pairs(amplicon_merged_seg,gene_properties,reciprocal_overlap = 0,svtype_matching = F)
map_amplicons_genes = map_amplicons_genes %>% dplyr::rename(patient_cna_merged_id=set1,gene_id=set2) 
map_amplicons_genes = map_amplicons_genes  %>% 
  left_join(amplicon_merged_segments_df[,c("patient_label","patient_cna_merged_id","amplicon_cr_l2fc_50","amplicon_footprint_id","amplicon_cna_width")],by="patient_cna_merged_id")
map_amplicons_genes = map_amplicons_genes %>% filter(overlap_set2_set1>0.33)

map_amplicons_genes = map_amplicons_genes %>% dplyr::mutate(amplicon_overlap = paste0(round(overlap_set2_set1,2), " at ",round(amplicon_cna_width/1e3),"kbp@",round(amplicon_cr_l2fc_50,2),"cr"))

# Patient gene level detailed  ----

  ## Load cohort gene cna
  #cohort_gene_cna = get_cohort_gene_cna(cohort,output_dir = utils_output_dir,baseline_shifts)
  cohort_gene_cna = get_gene_cna_from_segments(segments,segments_df,gene_properties,gene_properties_df)
  
  #if missing add tumors with genes
  cna_gene_lst = cohort_gene_cna[,c("gene_name","ensembl_id")] %>% unique()

  cohort_gene_cna_incl_missing_data = cohort_gene_cna
  
  all_genes_patients = tidyr::crossing(patient_label=cohort$patient_label, ensembl_id=cna_gene_lst$ensembl_id) 
  all_genes_patients %>% nrow() == (cna_gene_lst %>% nrow() * cohort %>% nrow())
  
  missing_genes=anti_join(all_genes_patients,cohort_gene_cna_incl_missing_data)
  
  ## NB: 2022-10-22 had to run this regardless of not matching
  #if((missing_genes %>% nrow() + cohort_gene_cna_incl_missing_data %>% nrow() ) == (cna_gene_lst %>% nrow() * cohort %>% nrow())) {
  print("Adding missing rows. check CN input data")
  missing_genes = missing_genes %>% 
    left_join(cna_gene_lst) %>%
    left_join(gene_properties_df[,names(gene_properties_df)[names(gene_properties_df)%in% names(cohort_gene_cna_incl_missing_data)]])
  
  cohort_gene_cna_incl_missing_data = rbind_no_colmatch(missing_genes,cohort_gene_cna_incl_missing_data)
  
  #}
  
  (cohort_gene_cna_incl_missing_data %>% nrow() ) == (cna_gene_lst %>% nrow() * cohort %>% nrow())
  cohort_gene_cna_incl_missing_data$patient_label %>% unique() %>% length() == cohort$patient_label %>% unique() %>% length()
  
  ### to here for cohort gene cna
  
  
## patient overview df for statistical tests ----
  #unfiltered_cohort_snv_somatic sv_gene_all_inside
  #not filtered by impact and whether they affect exons 
  

if(length(Sys.glob(patient_gene_overview_no_impact_filter_path))>0 & override_if_exists==F) {
  print("Output already exists, did not override")
  patient_gene_overview_no_impact = read.table(patient_gene_overview_no_impact_filter_path,sep="\t",header = T) 
} else {
    
  unfiltered_snvs = unfiltered_cohort_snv_somatic %>% dplyr::filter(gene_name!=""&!is.na(gene_name)) %>% 
    mutate(protein_consequence=ifelse(Consequence=="missense_variant",
                                      paste0("p.",str_replace(Amino_acids,"/",Protein_position)),
                                      Consequence)) %>%
    group_by(patient_label,gene_name) %>% 
    summarize(snv= lst_str(paste0(tumor_AF," taf (",protein_consequence," ",polyphen_label,"/",sift_label,")")),
              snv_cnt=cnt_str(patient_snv_id),
              snv_max_taf=max(tumor_AF))
  
  
  #Direct SNV SV bp alterations only
  patient_gene_overview_no_impact = cohort_gene_cna_incl_missing_data %>% 
    dplyr::select(patient_label,gene_name,ensembl_id,call,cr_l2fc_50,maf_50) %>% 
    left_join(unfiltered_snvs) %>%
    left_join(genes_sv_bp) 
  
  
  patient_gene_overview_no_impact = patient_gene_overview_no_impact %>%
    dplyr::mutate(across(.cols=contains("_cnt"),.fns = function(x){ifelse(!is.na(x),x,0)})) 
  
  patient_gene_overview_no_impact = patient_gene_overview_no_impact %>% filter(sv_bp_cnt>0|snv_cnt>0)
  
  patient_gene_overview_no_impact = patient_gene_overview_no_impact %>% 
    mutate(flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
           flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name) 
  
  write.table(patient_gene_overview_no_impact,patient_gene_overview_no_impact_filter_path,sep="\t",quote=F,row.names = F,col.names = T) 
}
  
  
## Patient detailed ----
if(length(Sys.glob(patient_gene_overview_detailed_path))>0 & override_if_exists==F) {
    print("Output already exists, did not override")
    patient_gene_overview_detailed = read.table(patient_gene_overview_detailed_path,sep="\t",header = T) 
} else {
    
  ## Detailed version will be wide dataframe for manual processing
  
  #checked, protein consequence also good for reverse strand genes
  ##2023-04-06 add hotspot promoter mutations => TERT
  snv_promoter_hotspot = unfiltered_cohort_snv_somatic %>% filter(snv_id %in% snv_promoter_hotspot_lst)
  
  cohort_snv_somatic_incl_promoter = rbind(cohort_snv_somatic,snv_promoter_hotspot)
  
  snvs = cohort_snv_somatic_incl_promoter %>% dplyr::filter(gene_name!=""&!is.na(gene_name)) %>% 
    mutate(protein_consequence=ifelse(Consequence=="missense_variant",
                                      paste0("p.",str_replace(Amino_acids,"/",Protein_position)),
                                      Consequence)) %>%
    group_by(patient_label,gene_name) %>% 
    summarize(snv= lst_str(paste0(tumor_AF," taf (",protein_consequence," ",polyphen_label,"/",sift_label,")")),
              snv_cnt=cnt_str(patient_snv_id),
              snv_max_taf=max(tumor_AF))

  #SV bp means sv breakpoint inside. $svtype $sv_tumor_af taf $sv_position $sv_length =>sv position is relative to the gene either inside or partially overlappin
  #NB: 
  # sv_gene_long == filtered_cohort_sv_gene_overlaps but summarized
  # svs_1mb_gene == flank_gene_overlaps anti joined with  filtered_cohort_sv_gene_overlaps to remove intersects so only 1mb nearby
  # 2023-08-22 update remove the patient sv merged not patient gene level 
  
  patient_gene_overview_detailed = cohort_gene_cna_incl_missing_data %>% 
    dplyr::select(patient_label,gene_name,ensembl_id,gene_id,call,cr_l2fc_50,maf_50) %>% 
    left_join(snvs) %>%
    left_join(genes_sv_bp) %>% 
    left_join(genes_sv_1mb) %>% 
    left_join(complex_genes_svs) %>% 
    left_join(complex_genes_svs_1mb) %>% 
    left_join(complex_genes_class)
  
  patient_gene_overview_detailed = patient_gene_overview_detailed %>%
    left_join(gene_properties_df[,c(gene_cols)] %>% unique()) %>%
    left_join(cohort[,c(cohort_anno_display_cols)]) 
  
  #relevant if  expression altered &&  sv spanning, 1mb sv, loss or gain 
  #always affected if snv, sv bp, > or expression altered and then look at these cols for the cause
  
  #need the or cohort level vst because not always CN altered 
  

  patient_gene_overview_detailed = patient_gene_overview_detailed %>%
    dplyr::mutate(across(.cols=contains("_cnt"),.fns = function(x){ifelse(!is.na(x),x,0)})) 
  
  ##amplification
  ## 2023-05-01 update: amplification relative to chromosome average 
   
  patient_gene_overview_detailed = patient_gene_overview_detailed %>% left_join(chrom_cn_fractions) %>% 
    dplyr::mutate(chrom_centered_cr_l2fc_50=(cr_l2fc_50-chrom_mean_cr))
  
 
  #add as amplicon overlap -> also flag_cn_amplification
  patient_gene_overview_detailed = patient_gene_overview_detailed %>% left_join(map_amplicons_genes[,c("gene_id","patient_label","amplicon_cr_l2fc_50","amplicon_overlap")]) %>%
    dplyr::mutate(flag_amplicon_overlap = !is.na(amplicon_overlap)) 
  
  ## 2023-07-05 update: amplification also via amplicon overlap
  patient_gene_overview_detailed = patient_gene_overview_detailed %>% 
    dplyr::mutate( flag_cn_amplification = flag_amplicon_overlap | 
                     (!is.na(chrom_centered_cr_l2fc_50) & (chrom_centered_cr_l2fc_50>cr_amplification_threshold & cr_l2fc_50>cr_amplification_threshold)))
  
  
  ## sv affected iff sv bp || (expression altered &&  sv spanning | 1mb sv)
  ## snv affected if snv
  ## cn affected if call!=0 && expression altered
  ## although high AF can also indicate both alleles hit > then often LOH observed 
  ## 2022-10-26 add amplification as type of affected by CNA 
  ## 2023-03-20 add homozygous loss as type of affected by CNA, no expression validation needed
  
  # todo later refactor flag_cn_hom_loss and flag_cn_amplification into call amp and hom_loss
  
  patient_gene_overview_detailed = patient_gene_overview_detailed %>% 
    mutate(
      flag_cn_hom_loss = !is.na(cr_l2fc_50) & cr_l2fc_50 < (cr_hom_loss_threshold)
    )
  
  
  ## later: maybe make it multi-modal hit? or just use it manually!
  # 2023-08-22 added flag_affected_by_sv_1mb_cn_change but not as "flag affected"
  patient_gene_overview_detailed = patient_gene_overview_detailed  %>%
    mutate(flag_affected_by_sv =  (sv_bp_cnt>0 ),
           flag_affected_by_cna = (flag_cn_amplification | flag_cn_hom_loss),
           flag_affected_by_snv =  (snv_cnt>0),
           flag_affected_by_sv_1mb_cn_change = ifelse(is.na(cr_l2fc_50),FALSE,
             (sv_1mb_cnt>0 & ( (chrom_centered_cr_l2fc_50>cr_minimum_cn_change & cr_l2fc_50>cr_minimum_cn_change) | (chrom_centered_cr_l2fc_50<(-cr_minimum_cn_change) & cr_l2fc_50<(-cr_minimum_cn_change))))),
           flag_affected = (flag_affected_by_cna|flag_affected_by_sv|flag_affected_by_snv))
  
  
  patient_gene_overview_detailed$alteration=""
  patient_gene_overview_detailed$alteration_simple=""
  
  patient_gene_overview_detailed = patient_gene_overview_detailed %>% 
    mutate(alteration = ifelse(snv_cnt>0,paste0(alteration, " snv"),alteration),
           alteration = ifelse(sv_bp_cnt>0,paste0(alteration, " sv_bp"),alteration),
           alteration = ifelse(flag_cn_amplification,paste0(alteration, " amp"),alteration), ## add special for amplification
           alteration = ifelse(flag_cn_hom_loss,paste0(alteration, " hom_loss"),alteration), ## add special for hom loss
           alteration = ifelse(flag_affected_by_cna & flag_cn_amplification==F & flag_cn_hom_loss==F,paste0(alteration, " ",call),alteration),
           alteration = ifelse(alteration==""&flag_affected_by_sv_1mb_cn_change,paste0(alteration, " sv_1mb_cna"),alteration)) %>% 
    
    mutate(alteration_simple = ifelse(snv_cnt>0,paste0(alteration_simple, " snv"),alteration_simple),
           alteration_simple = ifelse(flag_affected_by_sv,paste0(alteration_simple, " sv"),alteration_simple),
           alteration_simple = ifelse(flag_cn_amplification,paste0(alteration_simple, " amp"),alteration_simple),
           alteration_simple = ifelse(flag_cn_hom_loss,paste0(alteration_simple, " hom_loss"),alteration_simple),
           alteration_simple = ifelse(flag_affected_by_cna & flag_cn_amplification==F & flag_cn_hom_loss==F,paste0(alteration_simple, " cna"),alteration_simple),
           alteration_simple = ifelse(alteration_simple==""&flag_affected_by_sv_1mb_cn_change,paste0(alteration_simple, " sv_1mb_cna"),alteration_simple))
  
  patient_gene_overview_detailed[patient_gene_overview_detailed$alteration=="",c("alteration")]=NA
  patient_gene_overview_detailed[patient_gene_overview_detailed$alteration_simple=="",c("alteration_simple")]=NA
  
  
  patient_gene_overview_detailed = patient_gene_overview_detailed %>% 
    mutate(flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
           flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name) 
  }  
# export ----  
  write.table(patient_gene_overview_detailed,patient_gene_overview_detailed_path,sep="\t",quote=F,row.names = F,col.names = T) 
  
  patient_gene_overview_detailed$patient_label %>% unique() %>% length()
  
  #save a cancer genes only version for faster loading
  
  patient_gene_overview_detailed_cancer = patient_gene_overview_detailed %>% dplyr::filter(flag_gene_of_interest)
  
  write.table(patient_gene_overview_detailed_cancer,patient_gene_overview_detailed_cancer_path,sep="\t",quote=F,row.names = F,col.names = T) 
  
  #Direct SNV SV bp alterations only
  write.table(patient_gene_overview_detailed %>% filter(sv_bp_cnt>0|snv_cnt>0),patient_gene_overview_snv_sv_path,sep="\t",quote=F,row.names = F,col.names = T) 
  write.table(patient_gene_overview_detailed %>% filter(flag_affected),patient_gene_overview_affected_path,sep="\t",quote=F,row.names = F,col.names = T) 
  

  
  
  