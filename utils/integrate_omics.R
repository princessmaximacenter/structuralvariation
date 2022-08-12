

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

source(paste0(wdir,"wilms.conf"))

## supplementary / manuscript tables ----

cohort_sup_table_path="~/data/metadata/cohort.wilms_v3.supplementary_table_1.tsv"

chr_arm_sup_table_path = "~/results/wilms_v3/chr_arm_cna.wilms_v3.tsv"
mutation_burden_sup_table_path = "~/results/wilms_v3/mutation_burden.wilms_v3.tsv"
wilcox_significant_sup_table_path = "~/results/wilms_v3/wilcox_expression_cna_.wilms_v3.tsv"
metagenes_bp_sup_table_path = "~/results/wilms_v3/ex_cluster_biological_processes.wilms_v3.tsv"
recurrent_cn_sv_overexpressed_sup_table_path = "~/results/wilms_v3/recurrent_cn_sv_overexpressed.wilms_v3.tsv"
chr1q_overexpressed_genes_1q_tumors_cancer_sup_table_path = "~/results/wilms_v3/chr1q_gain_tumors_overexpressed_genes_1q.cancer.wilms_v3.tsv"
patient_gene_overview_detailed_cancer_sup_table_path = "~/results/wilms_v3/patient_gene_overview.cancer.wilms_v3.tsv"
wnt_pathway_ex1_sup_table_path = "~/results/wilms_v3/wnt_pathway_ex1.wilms_v3.tsv"
mutation_burden_statistics_sup_table_path="~/results/wilms_v3/mutation_burden.statistical_tests.wilms_v3.tsv"
somatic_svs_sup_table_path="~/results/wilms_v3/somatic_svs.wilms_v3.tsv"

## settings ----


sv_tumor_af_threshold=0.1
snv_tumor_af_threshold=0.1
autosomes= c(paste("chr",1:22,sep=""),"chrX")

l2fc_wilcox_threshold = log2(1.1)
fdr_wilcox_threshold = 0.2
vst_zscore_threshold=1.98
cohort_label="pretreated_bm_outlier"

utils_output_dir="~/results/utils/"

genes_of_interest_path = paste0(resources_dir,"genes_of_interest.tsv")

cohort_snv_path = "~/results/wilms_v2/cohort_snv.somatic.tsv"

cohort_gene_expression_path = "~/results/wilms_v3/cohort.gene_expression.tsv"
cohort_gene_cna_path = "~/results/wilms_v3/cohort.gene_cna.tsv"
cohort_gene_expression_zscore_path = paste0("~/results/wilms_v3/","gene_expression_zscore.",cohort_label,".tsv")
genes_vst_unfiltered_path= paste0("~/results/wilms_v3/unfiltered.",cohort_label,".blinded.","genes.vst.tsv")

patient_gene_overview_path = "~/results/wilms_v3/gene_level_overview.somatic.tsv"
patient_gene_overview_detailed_path = "~/results/wilms_v3/gene_level_overview.somatic.detailed.tsv"
patient_gene_overview_detailed_cancer_path = "~/results/wilms_v3/gene_level_overview.somatic.detailed.cancer.tsv"

gene_centric_results_dir= "~/results/wilms_v2/gene_centric/"
filtered_cohort_sv_gene_overlaps_path = paste0("~/results/wilms_v3/sv_gene_overlaps.somatic.tsv")
flank_gene_overlaps_path="~/results/wilms_v3/sv_gene.1mb_window.somatic.taf01.tsv"
somatic_svs_path = "~/results/wilms_v2/cohort_sv.somatic.anno.tsv"
merged_svs_path = "~/results/wilms_v3/merged_sv.somatic.tsv"
ctx_burden_path = "~/results/wilms_v3/ctx_burden.tsv"

cohort_chr_arm_cna_path = "~/results/wilms_v3/cohort_chr_arm_cna.tsv"
cohort_mutation_overview_path = "~/results/wilms_v3/cohort_mutation_burden.tsv"

expression_data_dir= "/Users/ianthevanbelzen/data/rna_expression/gencode31oct2019_multioverlap_largest_overlap/"
expression_gene_file_ext= "_RNA-Seq.gene_id.exon.counts.txt"

nmf_reports_dir="~/PycharmProjects/gene_expression/wilms/nmf/"
nmf_results_name="pretreated_bm_outlier_top10k_noX_blinded"
metagenes_path = paste0(nmf_reports_dir,nmf_results_name,"/metagenes.",nmf_results_name,".tsv" )


#expression + cn data, includes 1mb sv bp and neutral zscores too
## and adds wilcox flags
genes_cna_expression_path= "~/results/wilms_v3/cohort_gene_cna_expression.annotated.tsv"

wilcox_gain_path = paste0("~/results/wilms_v3/","expression_wilcox.","gain_vs_neutral",".",cohort_label,".tsv")
wilcox_loss_path = paste0("~/results/wilms_v3/","expression_wilcox.","loss_vs_neutral",".",cohort_label,".tsv")

wilcox_overall_path = paste0("~/results/wilms_v3/","expression_wilcox.","overall",".",cohort_label,".tsv")
wilcox_overall_wide_path = paste0("~/results/wilms_v3/","expression_wilcox.","overall_wide",".",cohort_label,".tsv")


cohort_anno_display_cols = c("patient_id","expression_cluster","cn_similarity_cluster","c1a_pattern","histology_consensus")


# Resources ----

## Cancer genes ----

cancer_genes = get_cancer_genes(resources_dir)
cancer_related=get_cancer_genes(resources_dir) %>% dplyr::rename(gene_name=gene_id)


wilms_genes_of_interest = read.table(genes_of_interest_path, sep="\t",header=T)
wilms_genes_of_interest_predisposition = wilms_genes_of_interest %>% filter(subset=="predisposition")


genes_of_interest = rbind(cancer_related[,c("gene_name")] %>% as.data.frame(),
                          wilms_genes_of_interest[,c("gene_name")] %>% as.data.frame())
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


## Genes ----

gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)

gene_properties_df = get_gene_properties_df(gtf_path)

gene_properties=GRanges(gene_properties_df)
names(gene_properties) = gene_properties$gene_id

genes_to_cytobands = get_genes_to_cytobands(gene_properties,chromosome_bands)

gene_properties_df = gene_properties_df %>% left_join(genes_to_cytobands)
gene_properties_df$gene_coord = gene_properties_df$to_coordinate
gene_properties_df$gene_name_display =  paste0(gene_properties_df$gene_name,"_",gene_properties_df$cytoband,"_",gene_properties_df$ensembl_id)
gene_cols=c("ensembl_id","gene_name","gene_type","gene_coord","seqnames","cytoband","gene_name_display")
gene_cols = gene_cols[gene_cols %in% names(gene_properties_df)]

#annotate
gene_properties_df = gene_properties_df %>% 
  mutate(flag_cancer_gene = gene_name %in% cancer_genes$gene_id,
         flag_wilms_gene_of_interest = gene_name %in% wilms_genes_of_interest$gene_name,
         flag_gene_of_interest = gene_name %in% genes_of_interest$gene_name) 

gene_properties_df = gene_properties_df %>%  left_join(cancer_genes[,c("gene_id","oncogene","tsg")],by=c("ensembl_id"="gene_id")) %>% dplyr::rename(flag_oncogene=oncogene,flag_tsg=tsg)

#check gene_properties_df[duplicated(gene_properties_df$gene_id),]

gene_properties_df_flags = names(gene_properties_df)[grepl("flag",names(gene_properties_df))]

for(flag in gene_properties_df_flags) {
  gene_properties_df[is.na(gene_properties_df[,flag]),flag]=F
}


genes_1mb = GRanges(gene_properties_df)
names(genes_1mb) =genes_1mb$gene_id
genes_1mb =  resize(genes_1mb, width = width(genes_1mb)+1e6, fix = "center")


# Cohort ----
##Load cohort and annotate 

cohort = read.table(cohort_sup_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]

cohort_for_expression = cohort  %>% filter(flag_use_expression)


# Read variants ----
## somatic SV gene overlaps ----

if(length(Sys.glob(filtered_cohort_sv_gene_overlaps_path))==1){
  filtered_cohort_sv_gene_overlaps=read.table(filtered_cohort_sv_gene_overlaps_path,header=T,sep="\t")
} else {
  cohort_sv_gene_overlaps = data.frame()
  for(pid in cohort$patient_id) {
    patient = cohort %>% filter(patient_id==pid)
    map_template_vars_patient = c('${merged_svs_dir}'=gene_centric_results_dir,'${patient_basename}'=patient$basename)  
    
    svs_union_gene_overlaps_anno_path = stri_replace_all_fixed(svs_union_gene_overlaps_anno_somatic_path_multitool_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    
    if(length(Sys.glob(svs_union_gene_overlaps_anno_path))==1) { 
      #Load sv gene overlaps gene centric
      sv_gene_overlaps = read.table(svs_union_gene_overlaps_anno_path,header=T,sep="\t",stringsAsFactors = F)
      cohort_sv_gene_overlaps = rbind(cohort_sv_gene_overlaps,sv_gene_overlaps)
    } else {
      print(paste0("Missing for patient ",patient$patient_id,": ",svs_union_gene_overlaps_anno_path))
      next()
    }  
    
  }        
  
  cohort_sv_gene_overlaps = cohort_sv_gene_overlaps %>% mutate(flag_cancer_gene = gene_name %in% cancer_genes$gene_id)
  cohort_sv_gene_overlaps = cohort_sv_gene_overlaps %>% mutate(flag_wilms_gene_of_interest = gene_name %in% wilms_genes_of_interest$gene_name )
  cohort_sv_gene_overlaps = cohort_sv_gene_overlaps  %>% mutate(flag_gene_of_interest = gene_name %in% genes_of_interest$gene_name )
  
  cohort_sv_gene_overlaps$tumor_af=as.numeric(cohort_sv_gene_overlaps$tumor_af)
  cohort_sv_gene_overlaps$overlap_gene_sv_frac=as.numeric(cohort_sv_gene_overlaps$overlap_gene_sv_frac)
  cohort_sv_gene_overlaps$overlap_sv_gene_frac=as.numeric(cohort_sv_gene_overlaps$overlap_sv_gene_frac)
  cohort_sv_gene_overlaps$gene_overlap_tool_cnt=as.numeric(cohort_sv_gene_overlaps$gene_overlap_tool_cnt)
  
  cohort_sv_gene_overlaps = cohort_sv_gene_overlaps %>% annotate_sv_affects_exon() %>% annotate_sv_gene_position() %>% annotate_sv_repeat_mediated() #%>% make_svlen_bin() %>% make_tumor_af_bin() 
  
  #filter by tumor_af>0.1 and exon is affected
  
  filtered_cohort_sv_gene_overlaps = cohort_sv_gene_overlaps %>% 
    filter(tumor_af>sv_tumor_af_threshold) %>%
    filter(tool_cnt>1 & gene_overlap_tool_cnt>1) %>% 
    filter(flag_affects_exon==T)
  
  write.table(filtered_cohort_sv_gene_overlaps, filtered_cohort_sv_gene_overlaps_path,sep="\t",quote = F,row.names = F,col.names = T)
}


sv_gene_long = filtered_cohort_sv_gene_overlaps %>%
  group_by(patient_id,gene_name,gene_id,sv_position,svtype,patient_sv_merged) %>% summarize(
    svlen_mean=mean(svLen,na.rm=T),
    tumor_af_mean=mean(tumor_af,na.rm=T)
  )


## somatic SNVs ----

##Flags cancer gene, wilms gene of interest > TODO flags from gene properties? 

cohort_snv_somatic = read.table(cohort_snv_path,sep="\t",header=T)
cohort_snv_somatic = cohort_snv_somatic %>% filter(patient_id %in% cohort$patient_id)

cohort_snv_somatic$patient_snv_id=paste0(cohort_snv_somatic$patient_id,"_",cohort_snv_somatic$snv_id)

cohort_snv_somatic = cohort_snv_somatic %>% mutate(dbsnp_present = grepl("rs",Existing_variation))
cohort_snv_somatic = cohort_snv_somatic %>% mutate(cosmic_present = grepl("COSV",Existing_variation))

cohort_snv_somatic = cohort_snv_somatic %>% mutate(flag_cancer_gene = gene_name %in% cancer_genes$gene_id )
cohort_snv_somatic = cohort_snv_somatic %>% mutate(flag_wilms_gene_of_interest = gene_name %in% wilms_genes_of_interest$gene_name )
cohort_snv_somatic = cohort_snv_somatic  %>% mutate(flag_gene_of_interest = gene_name %in% genes_of_interest$gene_name)

cohort_snv_somatic = cohort_snv_somatic %>% separate(col=PolyPhen,sep="\\(",into=c("polyphen_label","score"),remove=F)
cohort_snv_somatic = cohort_snv_somatic %>% separate(col=SIFT,sep="\\(",into=c("sift_label","score"),remove=F)

cohort_snv_somatic = cohort_snv_somatic %>% mutate(flag_protein_function_affected = ( grepl("damaging",PolyPhen) | grepl("deleterious",SIFT)))
cohort_snv_somatic = cohort_snv_somatic %>% mutate(flag_benign = ( grepl("benign",PolyPhen) | grepl("tolerated",SIFT)))

#always filter by tumor AF, autosomes+chrX, remove blacklisted , remove dbsnp if not in cosmic as well 
#cohort_snv_somatic_unfiltered = cohort_snv_somatic #cohort_snv_somatic_bk in other file

cohort_snv_somatic = cohort_snv_somatic %>% filter(tumor_AF>snv_tumor_af_threshold)
cohort_snv_somatic = cohort_snv_somatic %>% filter(seqnames %in% autosomes) %>% filter(!dac_blacklist)
#remove dbsnp unless also cosmic 
cohort_snv_somatic = cohort_snv_somatic %>% filter(!dbsnp_present | cosmic_present)

### DIVERGENCE OF GENE AFFECTING AND BURDEN
cohort_snv_somatic_first_filtering= cohort_snv_somatic 

# filtering for effect on gene 
#Filter by impact for if affects gene
cohort_snv_somatic = cohort_snv_somatic_first_filtering %>% filter(IMPACT %in% c("MODERATE","HIGH"))

#remove benign unless also damaging (conflict sift and polyphen) or cosmic
cohort_snv_somatic = cohort_snv_somatic %>% filter(!flag_benign | flag_protein_function_affected | cosmic_present)


## SNV burden ----
if(FALSE) {
denominator = 41

## checked what would be removed and decided to not filter burden by DP/AD only for the mutsig calling then
#filter_out_questionmark = cohort_snv_somatic %>% filter(tumor_DP<20|normal_DP<20 | !grepl(", 0)",normal_AD))
#filter_out_questionmark %>% filter(flag_gene_of_interest)

# NO additional filtering for burden
#filter DP tumor and normal >20, remove alt_supporting_read_in_normal
burden_cohort_snv_somatic = cohort_snv_somatic_first_filtering %>% filter(seqnames !="chrX")

## choose to remove LOW in its entirety because the splice affecting mutations are not in coding REGION 
burden_cohort_snv_somatic_nonsynonymous = burden_cohort_snv_somatic  %>% 
  filter(IMPACT %in% c("HIGH","MODERATE") & !Consequence %in% c("regulatory_region_ablation"))

burden_cohort_snv_somatic_nonsynonymous = burden_cohort_snv_somatic_nonsynonymous %>% left_join(genes_vep104_df[,c("ensembl_id","gene_biotype")])
burden_cohort_snv_somatic_nonsynonymous = burden_cohort_snv_somatic_nonsynonymous %>% filter( gene_biotype =="protein_coding")

burden_cohort_snv_somatic_nonsynonymous_only_snv =  burden_cohort_snv_somatic_nonsynonymous %>% filter(VARIANT_CLASS=="SNV")

cohort_snv_burden_somatic = burden_cohort_snv_somatic %>% cnt_rows(attr_name = "snv_indel_somatic_cnt") 
cohort_snv_burden_somatic = cohort_snv_burden_somatic %>% left_join( burden_cohort_snv_somatic_nonsynonymous %>% cnt_rows(attr_name = "snv_indel_somatic_nonsyn_cnt") ) 
cohort_snv_burden_somatic = cohort_snv_burden_somatic %>% left_join( burden_cohort_snv_somatic_nonsynonymous_only_snv %>% cnt_rows(attr_name = "snv_somatic_nonsyn_cnt") ) 
cohort_snv_burden_somatic[is.na(cohort_snv_burden_somatic)]=0


cohort_snv_burden_somatic$tmb = cohort_snv_burden_somatic$snv_indel_somatic_nonsyn_cnt/denominator

cohort_snv_burden_somatic$tmb %>% summary()
cohort_snv_burden_somatic$snv_indel_somatic_nonsyn_cnt %>% summary()
cohort_snv_burden_somatic$snv_somatic_nonsyn_cnt %>% summary()

#write.table(cohort_snv_burden_somatic,"~/results/wilms_v3/cohort_snv_burden_somatic.tsv",quote=F,sep="\t",row.names = F,col.names = T)
}
## somatic SVs ----

svs_df = read.table(somatic_svs_path,sep="\t",header=T)
svs_df$patient_sv_merged = paste0(svs_df$patient_id,"_",svs_df$sv_merged)
svs_df$partner_sv_name = paste0(svs_df$patient_id,"_",svs_df$partner) #for multitool ctx

svs_df = svs_df %>% filter(patient_id %in% cohort$patient_id) 

multitool_rs6 = get_multi_tool_support(svs_df)
svs_df = svs_df %>% filter(patient_sv_merged %in% multitool_rs6$patient_sv_merged)
svs_df = svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
svs_df = svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
unfiltered_svs_df = svs_df

svs_df = unfiltered_svs_df %>% filter(tumor_af>sv_tumor_af_threshold)

## for most coding analyses filter by tumor AF >0.1 and then reassess multi tool support 

filtered_svs_multi_tool_support = get_multi_tool_support(svs_df)
svs_df = svs_df %>% filter(patient_sv_merged %in% filtered_svs_multi_tool_support$patient_sv_merged) 

#filtered_out_svs = merged_svs_tumor_af %>% filter( tumor_af_mean>0.1 & !patient_sv_merged %in% svs_df$patient_sv_merged)
#filtered_out_svs %>% View()

svs = GRanges(svs_df)
names(svs) = svs$patient_sv_name

svs_df[duplicated(svs_df$patient_sv_name),] %>% nrow() ==0



unfiltered_svs_df = unfiltered_svs_df %>% mutate(flag_in_filtered_svs=patient_sv_name %in% svs_df$patient_sv_name)
#use flag to group because even merging can be different before/after filtering
merged_svs = unfiltered_svs_df %>%
  #   head() %>% 
  group_by(patient_id,svtype,patient_sv_merged,sv_merged_coordinate,flag_in_filtered_svs) %>% summarize(
    tumor_af_mean=mean(tumor_af,na.rm=T),
    svlen_mean=mean(svLen,na.rm=T),
    svlen_mean_str = ifelse(is.na(svlen_mean),"",
           ifelse(svlen_mean>1e6,paste0(round(svlen_mean/1e6,2)," Mbp"),
                  paste0(round(svlen_mean/1e3,2)," kbp"))),
    chrom=toString(sort(unique(seqnames))),
    partner_chrom=toString(sort(unique(partner_chrom))),
    chr_arm=toString(sort(unique(chr_arm))),
    partner_chr_arm=toString(sort(unique(partner_chr_arm))),
    cytoband=toString(sort(unique(cytoband))),
    partner_cytoband=toString(sort(unique(partner_cytoband))),
    coordinates=toString(coordinate),
    tools=toString(tool),
  )

merged_svs$cytoband = factor(merged_svs$cytoband,levels=gtools::mixedsort(unique(merged_svs$cytoband)))
merged_svs = merged_svs[order(merged_svs$cytoband),] %>% ungroup()

#export merged svs
merged_svs = merged_svs  %>% left_join(cohort[,cohort_anno_display_cols])  

#write.table(merged_svs,merged_svs_path,sep="\t",col.names = T,row.names = F,quote=F)

merged_svs %>% filter(grepl("chr1q",cytoband) & c1a_pattern !="no_c1a" & tumor_af_mean>0.1 & flag_in_filtered_svs==F) %>% select(patient_id,svtype,svlen_mean_str,tumor_af_mean)# %>% View()

#write.table(merged_svs %>% filter(patient_id %in% rloops_patients),"/Users/ianthevanbelzen/surfdrive/Shared/Rloops and SVs/merged_svs.somatic.anno.rloops.tsv",sep="\t",col.names = T,row.names = F,quote=F)
#write.table(svs_df %>% filter(patient_id %in% rloops_patients),"/Users/ianthevanbelzen/surfdrive/Shared/Rloops and SVs/svs.somatic.anno.rloops.tsv",sep="\t",col.names = T,row.names = F,quote=F)


## CTX burden  ----

svs_ctx = svs_df %>% filter(svtype=="CTX")

svs_ctx = svs_ctx %>% left_join(cohort_anno[,c("patient_id","cn_similarity_cluster","c1a_pattern")])  %>% unique()


ctx_burden = svs_ctx[,c("patient_id","partner_chrom","seqnames","patient_sv_merged")] %>% 
  filter(seqnames %in% autosomes & partner_chrom %in% autosomes) %>% 
  unique() %>% 
  group_by(patient_id,seqnames,partner_chrom) %>% summarize(ctx_per_chrom_pair=length(unique(patient_sv_merged)))

ctx_burden_svs = svs_ctx %>% 
  filter(seqnames %in% autosomes & partner_chrom %in% autosomes) %>% 
  group_by(patient_id) %>% summarize(ctx_svs_cnt=length(unique(patient_sv_merged)))

ctx_burden = ctx_burden %>%  group_by(patient_id) %>%
  summarize(ctx_chrom_cnt=length(unique(c(seqnames, partner_chrom))),
            links=toString(unique(sort( paste0(seqnames,"-",partner_chrom," (",ctx_per_chrom_pair,")")))),
            links_short=str_replace_all(toString(unique(sort( paste0("t(",seqnames,";",partner_chrom,")")))),"chr",""))

ctx_burden = ctx_burden_svs %>% left_join(ctx_burden) %>% add_missing_patients_to_overview(cohort)

ctx_burden[is.na(ctx_burden$ctx_svs_cnt),c("ctx_svs_cnt")]=0
ctx_burden[is.na(ctx_burden$ctx_chrom_cnt),c("ctx_chrom_cnt")]=0

#write.table(ctx_burden,ctx_burden_path,quote=F,sep="\t",row.names = F,col.names = T)
#ctx_burden = read.table(ctx_burden_path,sep="\t",header = T)




## SV burden  ----

source(paste0(utils_script_dir,"functions.sv_burden.R"))

group_cols=c("patient_id","svtype")
sv_merge_cols=c("patient_id","patient_sv_merged","svtype","sv_merged_coordinate") 

unfiltered_svs_merged_df = unfiltered_svs_df %>% filter(seqnames %in% autosomes & seqnames!="chrX") %>% dplyr::select(all_of(c("tool",sv_merge_cols))) %>% unique()
unfiltered_merged_sv_cnt_long = get_sv_burden_long(unfiltered_svs_merged_df,group_cols=group_cols,cnt_type="merged_row")
unfiltered_merged_sv_event_long = unfiltered_merged_sv_cnt_long %>% convert_cnt_sv_row_event()
cohort_sv_event_burden_somatic = unfiltered_merged_sv_event_long %>% filter(cnt_type=="merged_event" & tool=="all_cnt" & svtype=="all")

#filter by tumor af >0.1 like the snvs
filtered_svs_merged_df = svs_df %>% filter(seqnames %in% autosomes & seqnames!="chrX") %>% dplyr::select(all_of(c("tool",sv_merge_cols))) %>% unique()
merged_sv_cnt_long_taf = get_sv_burden_long(filtered_svs_merged_df,group_cols=group_cols,cnt_type="merged_row")
merged_sv_event_long_taf  = merged_sv_cnt_long_taf %>% convert_cnt_sv_row_event()
cohort_sv_event_burden_somatic_taf = merged_sv_event_long_taf %>% filter(cnt_type=="merged_event" & tool=="all_cnt" & svtype=="all")


# note 2022-04-25 in wilms v2 because all 32 patients so not the v3 cohort
# note 2022-06-19  v3 cohort excl chr X to make consistent with SNV burden
unfiltered_merged_sv_event_long$filter=NA
merged_sv_event_long_taf$filter="taf>0.1"

unfiltered_merged_sv_event_long$variable_name="somatic_sv_cnt"
merged_sv_event_long_taf$variable_name="somatic_sv_taf_cnt"

export_sv_event_long = rbind(unfiltered_merged_sv_event_long,merged_sv_event_long_taf)
#write.table(export_sv_event_long,"~/results/wilms_v3/sv_burden_merged_event.tsv",sep="\t",quote=F,row.names = F,col.names = T)
#export_sv_event_long = read.table("~/results/wilms_v3/sv_burden_merged_event.tsv",sep="\t",header = T)

#export_sv_event_long %>% filter(cnt_type=="merged_event" & tool=="all_cnt" & svtype=="all") %>%
#  select(patient_id,variable_name,cnt) %>% pivot_wider(names_from="variable_name",values_from="cnt")






## SVs 1Mb of genes flank_gene_overlaps ----


if(length(Sys.glob(flank_gene_overlaps_path))==1){
  flank_gene_overlaps=read.table(flank_gene_overlaps_path,header=T,sep="\t")
} else {
  
  flank_gene_overlaps = get_sv_gene_flank_overlaps(svs,genes_1mb,svs_df,gene_properties_df)
  
  flank_gene_overlaps = flank_gene_overlaps %>% annotate_gene_group(cancer_genes,wilms_genes_of_interest)
  
  flank_gene_overlaps = flank_gene_overlaps %>% separate(col=patient_sv_name,sep="_",into=c("patient_id","sv_name"),remove=F)
  flank_gene_overlaps = flank_gene_overlaps %>% filter(patient_id %in% cohort$patient_id)
  
  #dont use this function because of ctx multitool_support_svs_1mb_gene = get_multi_tool_support(svs_1mb_gene)
  multitool_support_svs_1mb_gene = flank_gene_overlaps %>% filter(grepl("merged",patient_sv_merged)) %>% 
    group_by(patient_sv_merged,gene_id) %>% summarize(gene_1mb_tool_cnt=length(unique(tool)))
  multitool_support_svs_1mb_gene = filter(multitool_support_svs_1mb_gene,gene_1mb_tool_cnt>1)
  
  flank_gene_overlaps = flank_gene_overlaps %>% left_join(multitool_support_svs_1mb_gene)
  
  write.table(flank_gene_overlaps,flank_gene_overlaps_path,quote=F,sep="\t",row.names = F,col.names = T)
  
}

#Reapply multi tool filtering require same gene by tool
flank_gene_overlaps = flank_gene_overlaps %>% filter(gene_1mb_tool_cnt>1)

#remove sv bp intersects for gene centric analysis > not for expression analysis
#cannot do: #mutate(parsed_as_sv_bp_in_gene = patient_sv_merged %in% filtered_cohort_sv_gene_overlaps$patient_sv_merged) %>% 
#because then it doesnt take into account gene match

svs_1mb_gene = flank_gene_overlaps %>% 
  anti_join(filtered_cohort_sv_gene_overlaps[,c("patient_sv_merged","gene_id")],by=c("patient_sv_merged","gene_id"))

summary_svs_1mb_gene = svs_1mb_gene %>% 
  group_by(patient_id,gene_name,svtype,patient_sv_merged) %>% summarize(
    svlen_mean=mean(svLen,na.rm=T),
    tumor_af_mean=mean(tumor_af,na.rm=T),
    distance_kbp_mean=mean(distance_kbp,na.rm=T)
  )



## gene CNA ----
if(length(Sys.glob(cohort_gene_cna_path))==1) {
  cohort_gene_cna = read.table(cohort_gene_cna_path,sep="\t",header = T)
} else{
  
  cohort_gene_cna= data.frame()
  
  for(pid in cohort$patient_id) {
    patient = cohort %>% filter(patient_id == pid)
    
    gene_cna_path = paste0(utils_output_dir,gene_cna_outfile,patient$tumor_id,".tsv")
    gene_cna =read.table(gene_cna_path,sep="\t",header=T)
    cohort_gene_cna =rbind(cohort_gene_cna,gene_cna)
  }
  
  write.table(cohort_gene_cna,cohort_gene_cna_path,sep="\t",row.names = F,col.names = T,quote=F)
}

cohort_gene_cna = cohort_gene_cna %>% filter(patient_id %in% cohort$patient_id)
cohort_gene_cna$patient_id %>% unique() %>% length()
cohort_gene_cna$ensembl_id = remove_version_from_id(cohort_gene_cna$gene_id)
cohort_gene_cna = cohort_gene_cna %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4)

# Read expression ----

if(length(Sys.glob(cohort_gene_expression_zscore_path))==1) {
  cohort_gene_expression = read.table(cohort_gene_expression_zscore_path,sep="\t",header = T)
} else{
  cohort_gene_expression = data.frame(stringsAsFactors = F)
  for(pid in cohort_for_expression$patient_id) {
    patient = cohort %>% filter(patient_id == pid)
    gene_expr_filepath = paste0(expression_data_dir,patient$rna_id,expression_gene_file_ext)
    gene_expression_df = read_expression_data(gene_expr_filepath)
    gene_expression_df$patient_id = patient$patient_id
    
    cohort_gene_expression = rbind(cohort_gene_expression,gene_expression_df[,c("ensembl_id","gene_id","gene_name","counts","cpm","fpkm","patient_id")])
  }
  
  cohort_gene_expression$ensembl_id=factor(cohort_gene_expression$ensembl_id)
  cohort_gene_expression_summary = cohort_gene_expression %>% group_by(ensembl_id) %>% 
    mutate(sum_counts= sum(counts))
  cohort_gene_expression_summary[cohort_gene_expression_summary$sum_counts>10,] %>% nrow()
  cohort_gene_expression_summary[cohort_gene_expression_summary$sum_counts<10,] %>% nrow()
  
  #remove lowly expressed across samples  
  cohort_gene_expression = cohort_gene_expression %>% filter(ensembl_id %in% filter(cohort_gene_expression_summary,sum_counts>10)$ensembl_id)
  
  
  cohort_gene_expression = calculate_fpkm_zscore(cohort_gene_expression)
  
  write.table(cohort_gene_expression,cohort_gene_expression_zscore_path,sep="\t",row.names = F,col.names = T,quote=F)
}


# Read metagenes ----
metagenes = read.table(metagenes_path,sep="\t",header=T)

ex_clusters = metagenes$expression_cluster %>% unique()
metagenes_top_x=data.frame()
top_n=50
for(ex_cluster in ex_clusters) {
  top_metagenes = metagenes %>% filter(expression_cluster==ex_cluster) %>% arrange(-abs(l2fc)) %>% head(n=top_n)
  metagenes_top_x = rbind(metagenes_top_x,top_metagenes)
}

map_metagenes_top_x = metagenes_top_x %>% group_by(ensembl_id) %>% summarize(expr_cluster = toString(sort(unique(expression_cluster))))
metagenes$expr_cluster=metagenes$expression_cluster

metagenes = metagenes %>% mutate(flag_top_50 = ensembl_id %in% map_metagenes_top_x$ensembl_id)


# Read/generate normalized counts for expression analysis ----

if(length(Sys.glob(genes_vst_unfiltered_path))!=1) {
  
  #generate vst without the rowmeans >10 filter such that more genes are included. > this is different than for the clustering but am not going to repeat that now
  # blinded design 
  library(DESeq2)
  
  expression_df = read.table(cohort_gene_expression_path,sep="\t",header=T)
  expression_overview_wide = prepare_expression_df_deseq2(expression_df,as.character(cohort$patient_id),filter_read_counts = 0)
  
  if(length(names(expression_overview_wide)) !=  length(cohort$patient_id)) {
    print("Mismatch expression overview and cohort")
    #quit()
  }
  
  ## Make DDS object ----
  annotation=cohort
  rownames(annotation) = cohort$patient_id
  dds = DESeqDataSetFromMatrix(countData = expression_overview_wide,
                               colData = annotation,
                               design = ~ 1 )
  
  #For wilms cohort I used at least 5 in 3 samples
  nrow(dds)
  keep <- rowSums(counts(dds) >= 5) >= 3
  dds <- dds[keep,]
  nrow(dds)
  
  dds =DESeq(dds)
  
  #blinded
  vst = varianceStabilizingTransformation(dds,blind=T)
  
  fullset_genes = as.data.frame(assay(vst))
  fullset_genes$ensembl_id = rownames(fullset_genes)
  write.table(fullset_genes,genes_vst_unfiltered_path,quote=F,sep="\t",row.names = F)
} else {
  load_vst(genes_vst_unfiltered_path,as_global=T)
}


# CN vs expression ----

## Make df with expression and CN  ----
cohort_gene_cna_for_expression =  cohort_gene_cna %>% filter(patient_id %in% cohort_for_expression$patient_id)

genes_cna_expression = cohort_gene_cna_for_expression %>% left_join(cohort_gene_expression) 
genes_cna_expression = genes_cna_expression %>% left_join(gene_properties_df[,c("gene_id",gene_cols,gene_properties_df_flags)])

for(flag in gene_properties_df_flags) {
  genes_cna_expression[is.na(genes_cna_expression[,flag]),flag]=F
}


incomplete_genes_lst = genes_cna_expression[is.na(genes_cna_expression$fpkm),c("gene_id")]

genes_cna_expression = genes_cna_expression %>% filter(!gene_id %in% incomplete_genes_lst) 

#Annotate with cohort cols
genes_cna_expression = genes_cna_expression %>% left_join(cohort[,c("patient_id","c1a_pattern","histology_consensus","expression_cluster")])
genes_cna_expression = genes_cna_expression %>% mutate(flag_1q_gain = c1a_pattern != "no_c1a")

#annotate with metagenes and 10 k most var 
genes_cna_expression = genes_cna_expression  %>% left_join(metagenes[,c("ensembl_id","expr_cluster","flag_top_50")])
genes_cna_expression[is.na(genes_cna_expression$flag_top_50),c("flag_top_50")]=F

genes_cna_expression = genes_cna_expression %>% mutate(most_var_10k = !is.na(expr_cluster))


#add vst normalized counts to gene expression dataframe
genes_vst_long =  genes_vst %>% as.data.frame()
genes_vst_long$ensembl_id = rownames(genes_vst_long)
genes_vst_long =  genes_vst_long %>% pivot_longer(-ensembl_id,names_to = "patient_id",values_to="vst")
genes_cna_expression = genes_cna_expression %>% left_join(genes_vst_long)



#add info on SVs to genes_cna_expression
#flag sv bp > all sv nearby and inside
flank_gene_overlaps$flag_sv_bp=T
genes_cna_expression = genes_cna_expression %>% left_join(unique(flank_gene_overlaps[,c("patient_id","gene_id","flag_sv_bp")]),by=c("patient_id","gene_id"))
genes_cna_expression[is.na(genes_cna_expression$flag_sv_bp),c("flag_sv_bp")]=F

#joined gain/loss and 1mb sv
genes_cna_expression = genes_cna_expression %>% mutate(
  flag_sv_bp_or_gain  = ifelse((flag_sv_bp==T|call=="gain")&call!="loss","sv_1mb_or_gain",call),
  flag_sv_bp_or_loss  = ifelse((flag_sv_bp==T|call=="loss")&call!="gain","sv_1mb_or_loss",call))

#add overlapping sv / spanning and partial and inside
sv_gene_long$flag_sv_overlap=T
genes_cna_expression = genes_cna_expression %>% left_join(unique(sv_gene_long[,c("patient_id","gene_id","flag_sv_overlap")]),by=c("patient_id","gene_id"))
genes_cna_expression[is.na(genes_cna_expression$flag_sv_overlap),c("flag_sv_overlap")]=F

#joined gain/loss and sv overlap
genes_cna_expression = genes_cna_expression %>% mutate(
  flag_sv_and_gain  = ifelse((flag_sv_overlap==T|flag_sv_bp==T)&call=="gain","sv_and_gain",call),
  flag_sv_and_loss  = ifelse((flag_sv_overlap==T|flag_sv_bp==T)&call=="loss","sv_and_loss",call))


## Calculate neutral zcores ----
# goal; be able to compare patient with loss/gain/mutation against neutral-only averages 


## Neutal z scores with vst counts 
neutral_zscores_vst = genes_cna_expression %>% filter(!is.na(vst)) %>% 
  filter(call=="0") %>% dplyr::select(ensembl_id,vst,gene_name) %>% unique() %>%
  group_by(ensembl_id,gene_name) %>% 
  summarize(neutral_vst_mean = mean(vst,na.rm = T), 
            neutral_vst_sd = sd(vst,na.rm = T))

genes_cna_expression = genes_cna_expression %>% left_join(neutral_zscores_vst)

genes_cna_expression = genes_cna_expression %>% mutate(neutral_vst_zscore = (vst - neutral_vst_mean )/ neutral_vst_sd)

#genes_cna_expression = genes_cna_expression %>% dplyr::select(-vst,-neutral_vst_zscore,-neutral_vst_mean,-neutral_vst_sd)

#regular vst zscore
genes_cna_expression_bk = genes_cna_expression
genes_cna_expression = genes_cna_expression %>% calculate_vst_zscore()

##  saved later with wilcox added and chr arms
#write.table(genes_cna_expression,genes_cna_expression_path,sep="\t",col.names = T,row.names = F)


## Wilcoxon test for gain/loss/1mb sv vs neutral ----

### 2022-06-22 redo the wilcox  based on vst normalized counts => archived the old files and use same paths as before
#implemented BH padj 

excl_gene_lst = filter(genes_cna_expression,is.na(vst))$gene_id %>% unique()
assessed_genes_expression_lst = filter(genes_cna_expression,!is.na(vst))$gene_id %>% unique()

#prevent running accidentally 

if(FALSE) {
  
#calculate recurrence => 2022-06-22 update - only for expressed genes and use threshold of 3 or more patients
recurrence_patient_threshold=2
recurrent_gene_cna = genes_cna_expression %>% filter(call!="0") %>% filter(gene_id %in% assessed_genes_expression_lst) %>% get_recurrent_gene_cna() %>% filter(patient_cnt>recurrence_patient_threshold)

## Gain
recurrently_gained = recurrent_gene_cna %>% filter(call=="gain")

expression_wilcoxon_gain = get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrently_gained$gene_id,gene_id_col = "gene_id", 
                                                            attr_col="call", attr_value="gain",attr_value_contrast="0",prefix="gain_")
expression_wilcoxon_gain = expression_wilcoxon_gain %>%  left_join(gene_properties_df[,c("gene_id",gene_cols)])
expression_wilcoxon_gain$BH= p.adjust(expression_wilcoxon_gain$pvalue, method = 'BH')
write.table(expression_wilcoxon_gain,wilcox_gain_path,sep="\t",col.names = T,row.names = F,quote = F)

## Loss
recurrently_lost = recurrent_gene_cna %>% filter(call=="loss") 

expression_wilcoxon_loss= get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrently_lost$gene_id,gene_id_col = "gene_id", 
                                                        attr_col="call", attr_value="loss",attr_value_contrast="0",prefix="loss_")

expression_wilcoxon_loss = expression_wilcoxon_loss %>%  left_join(gene_properties_df[,c("gene_id",gene_cols)])
expression_wilcoxon_loss$BH= p.adjust(expression_wilcoxon_loss$pvalue, method = 'BH')
write.table(expression_wilcoxon_loss,wilcox_loss_path,sep="\t",col.names = T,row.names = F,quote = F)

## SV bp flanking 1 Mb vs  expression 
# filtered by same tool and >0.1 taf
#checked sv bp anno and we're all good
if(FALSE){
  id_cols=c("patient_id","gene_id")
  attr_test_col="flag_sv_bp"
  
  genes_cna_expression[,c(id_cols,attr_test_col)] %>% filter(flag_sv_bp) %>% anti_join(flank_gene_overlaps[,c(id_cols,attr_test_col)])
  flank_gene_overlaps[,c(id_cols,attr_test_col)] %>% merge(genes_cna_expression[,c(id_cols)]) %>% anti_join(filter(genes_cna_expression,flag_sv_bp)[,c(id_cols,attr_test_col)])
  
}

#Assess only those in genes_cna_expression
recurrent_gene_1mb = genes_cna_expression %>% filter(flag_sv_bp==T) %>% filter(gene_id %in% assessed_genes_expression_lst) %>% cnt_patients()
recurrently_1mb_sv = recurrent_gene_1mb %>% filter(patient_cnt>recurrence_patient_threshold) 

expression_wilcoxon_sv_1mb = get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrently_1mb_sv$gene_id,gene_id_col = "gene_id", 
                             attr_col="flag_sv_bp", attr_value="TRUE",attr_value_contrast="FALSE",prefix="sv_1mb_")

expression_wilcoxon_sv_1mb = expression_wilcoxon_sv_1mb %>% left_join(gene_properties_df[,c("gene_id",gene_cols)])
expression_wilcoxon_sv_1mb$BH= p.adjust(expression_wilcoxon_sv_1mb$pvalue, method = 'BH')
write.table(expression_wilcoxon_sv_1mb,wilcox_sv_1mb_path,sep="\t",col.names = T,row.names = F,quote = F)


recurrent_gene_sv_1mb_or_gain = genes_cna_expression %>% filter(gene_id %in% assessed_genes_expression_lst) %>% 
  filter(flag_sv_bp_or_gain=="sv_1mb_or_gain") %>% cnt_patients()
recurrent_gene_sv_1mb_or_gain = recurrent_gene_sv_1mb_or_gain %>% filter(patient_cnt>recurrence_patient_threshold)


expression_wilcoxon_sv_1mb_or_gain = get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrent_gene_sv_1mb_or_gain$gene_id,gene_id_col = "gene_id", 
                                                          attr_col="flag_sv_bp_or_gain", attr_value="sv_1mb_or_gain",attr_value_contrast="0",prefix="sv_1mb_or_gain_")
expression_wilcoxon_sv_1mb_or_gain = expression_wilcoxon_sv_1mb_or_gain %>% left_join(gene_properties_df[,c("gene_id",gene_cols)])
expression_wilcoxon_sv_1mb_or_gain$BH= p.adjust(expression_wilcoxon_sv_1mb_or_gain$pvalue, method = 'BH')
write.table(expression_wilcoxon_sv_1mb_or_gain,wilcox_sv_1mb_or_gain_path,sep="\t",col.names = T,row.names = F,quote = F)



recurrent_gene_sv_1mb_or_loss = genes_cna_expression %>% filter(gene_id %in% assessed_genes_expression_lst) %>% 
  filter(flag_sv_bp_or_loss=="sv_1mb_or_loss") %>% cnt_patients()
recurrent_gene_sv_1mb_or_loss = recurrent_gene_sv_1mb_or_loss %>% filter(patient_cnt>recurrence_patient_threshold)


expression_wilcoxon_sv_1mb_or_loss = get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrent_gene_sv_1mb_or_loss$gene_id,gene_id_col = "gene_id", 
                                                               attr_col="flag_sv_bp_or_loss", attr_value="sv_1mb_or_loss",attr_value_contrast="0",prefix="sv_1mb_or_loss_")

expression_wilcoxon_sv_1mb_or_loss = expression_wilcoxon_sv_1mb_or_loss %>% left_join(gene_properties_df[,c("gene_id",gene_cols)])
expression_wilcoxon_sv_1mb_or_loss$BH= p.adjust(expression_wilcoxon_sv_1mb_or_loss$pvalue, method = 'BH')
write.table(expression_wilcoxon_sv_1mb_or_loss,wilcox_sv_1mb_or_loss_path,sep="\t",col.names = T,row.names = F,quote = F)

}
## SV overlap vs expression
# filtered by same tool and >0.1 taf
if(FALSE){
  #strict, not nearby but overlap => not used
  recurrence_patient_threshold=2
  
  recurrent_gene_sv_and_gain = genes_cna_expression %>% filter(gene_id %in% assessed_genes_expression_lst) %>% 
    filter(flag_sv_and_gain=="sv_and_gain") %>% cnt_patients()
  
  recurrent_gene_sv_and_gain = recurrent_gene_sv_and_gain %>% filter(patient_cnt>recurrence_patient_threshold)
  recurrent_gene_sv_and_gain %>% left_join(gene_properties_df[,c("gene_id",gene_cols)]) %>% group_by(seqnames) %>% count()
  
  expression_wilcoxon_sv_and_gain = get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrent_gene_sv_and_gain$gene_id,gene_id_col = "gene_id", 
                                                                    attr_col="flag_sv_and_gain", attr_value="sv_and_gain",attr_value_contrast="0",prefix="sv_and_gain_")
  
  expression_wilcoxon_sv_and_gain = expression_wilcoxon_sv_and_gain %>% left_join(gene_properties_df[,c("gene_id",gene_cols)])
  expression_wilcoxon_sv_and_gain$BH= p.adjust(expression_wilcoxon_sv_and_gain$pvalue, method = 'BH')
  
  write.table(expression_wilcoxon_sv_and_gain,wilcox_sv_and_gain_path,sep="\t",col.names = T,row.names = F,quote = F)
  
  recurrent_gene_sv_and_loss = genes_cna_expression %>% filter(gene_id %in% assessed_genes_expression_lst) %>% 
    filter(flag_sv_and_loss=="sv_and_loss") %>% cnt_patients()
  
  recurrent_gene_sv_and_loss = recurrent_gene_sv_and_loss %>% filter(patient_cnt>recurrence_patient_threshold)
  #recurrent_gene_sv_and_loss %>% left_join(gene_properties_df[,c("gene_id",gene_cols)]) %>% group_by(seqnames) %>% count()
  
  expression_wilcoxon_sv_and_loss = get_wilcox_expression_groups_vst(genes_cna_expression,target_gene_lst = recurrent_gene_sv_and_loss$gene_id,gene_id_col = "gene_id", 
                                                                 attr_col="flag_sv_and_loss", attr_value="sv_and_loss",attr_value_contrast="0",prefix="sv_and_loss_")
  
  expression_wilcoxon_sv_and_loss = expression_wilcoxon_sv_and_loss %>% left_join(gene_properties_df[,c("gene_id",gene_cols)])
  expression_wilcoxon_sv_and_loss$BH= p.adjust(expression_wilcoxon_sv_and_loss$pvalue, method = 'BH')
  
  write.table(expression_wilcoxon_sv_and_loss,wilcox_sv_and_loss_path,sep="\t",col.names = T,row.names = F,quote = F)
  
}

# Collect wilcox expression & parse for interpretation ----

## also annotates and makes new export of 
#genes_cna_expression

if(length(Sys.glob(wilcox_overall_path))==1) {
    
  wilcox_significant_all = read.table(wilcox_overall_path,sep="\t",header = T)
  wilcox_significant_wide = read.table(wilcox_overall_wide_path,sep="\t",header = T)

} else {

#Load wilcox results and filter default to FDR <0.2 fold change 1.1 median or mean log2(1.1) as threshold

expression_wilcoxon_loss = read.table(wilcox_loss_path,sep="\t",header=T)
expression_wilcoxon_gain = read.table(wilcox_gain_path,sep="\t",header=T)


expression_wilcoxon_gain = expression_wilcoxon_gain %>%
  filter(BH<fdr_wilcox_threshold) %>% 
  filter(l2fc_median >l2fc_wilcox_threshold | l2fc_mean > l2fc_wilcox_threshold) %>%
  arrange(-l2fc_mean) 

expression_wilcoxon_loss = expression_wilcoxon_loss %>% 
  filter(BH<fdr_wilcox_threshold) %>% 
  filter(l2fc_median <(-l2fc_wilcox_threshold) | l2fc_mean < (-l2fc_wilcox_threshold)) %>%
  arrange(-l2fc_mean) 

#Collect all wilcox results in single df 

expression_wilcoxon_gain$analysis="gain_vs_neutral"
expression_wilcoxon_loss$analysis="loss_vs_neutral"


wilcox_significant_all = rbind_no_colmatch(expression_wilcoxon_gain,expression_wilcoxon_loss)

wilcox_significant_wide = wilcox_significant_all %>% select("ensembl_id","gene_id","gene_name_display","analysis") %>% mutate(flag=T)
wilcox_significant_wide = wilcox_significant_wide %>% pivot_wider(names_from="analysis",values_from="flag")
wilcox_significant_wide[is.na(wilcox_significant_wide)]=F

wilcox_significant_wide = wilcox_significant_wide %>% dplyr::rename_with(function(x){paste0("wilcox_",x)},!contains("gene_")&!contains("ensembl_id"))

wilcox_significant_wide$flag_any_wilcox=T


##annotate and export wilcox ----

wilcox_significant_wide = wilcox_significant_wide %>%  left_join(gene_properties_df[gene_properties_df$seqnames!="chrY",c(gene_cols,gene_properties_df_flags)])
wilcox_significant_all = wilcox_significant_all %>% left_join(gene_properties_df[gene_properties_df$seqnames!="chrY",c(gene_cols,gene_properties_df_flags)])


for(flag in c(gene_properties_df_flags)) {
  wilcox_significant_all[is.na(wilcox_significant_all[,flag]),flag]=F
  wilcox_significant_wide[is.na(wilcox_significant_wide[,flag]),flag]=F
}


write.table(wilcox_significant_all,wilcox_overall_path,sep="\t",col.names = T,row.names = F,quote = F)
write.table(wilcox_significant_wide,wilcox_overall_wide_path,sep="\t",col.names = T,row.names = F,quote = F)

}

## Add wilcox resuts to genes_cna_expression dataframe ----
#prevent auto rerun
if(FALSE){
genes_cna_expression = genes_cna_expression %>% left_join(wilcox_significant_wide)
genes_cna_expression = genes_cna_expression %>% dplyr::mutate(across(contains("wilcox"),~ replace_na(.,F)))

#only keep flags if patient has corresponding alteration
genes_cna_expression = genes_cna_expression %>% mutate(wilcox_gain_vs_neutral = ifelse(call=="gain",wilcox_gain_vs_neutral,F),
                                                       wilcox_loss_vs_neutral = ifelse(call=="loss",wilcox_loss_vs_neutral,F))
genes_cna_expression = genes_cna_expression %>% mutate(flag_any_wilcox=if_any(contains("wilcox_")))
#write.table(genes_cna_expression,genes_cna_expression_path,sep="\t",col.names = T,row.names = F)

}
#%>% filter(if_any(.cols=contains("wilcox"))


#  Chromosomes ----

if(length(Sys.glob(cohort_chr_arm_cna_path))==1) { 
  cohort_chr_arm_cna = read.table(cohort_chr_arm_cna_path,sep="\t",header = T)
  
} else {
cohort_chr_arm_cna= data.frame()
#cohort_fga_per_chrom = data.frame()
utils_output_dir="~/results/utils/"

for(pid in cohort$patient_id) {
  patient = cohort %>% filter(patient_id == pid)
  map_template_vars = c('${utils_output_dir}'=utils_output_dir,'${biomaterial_id}'=patient$tumor_id)
  chr_arm_cna_path = stri_replace_all_fixed(chr_arm_cna_path_template,names(map_template_vars), map_template_vars,vectorize=F)
  fga_per_chrom_path = stri_replace_all_fixed(fga_per_chrom_path_template,names(map_template_vars), map_template_vars,vectorize=F)
  
  chr_arm_cna =read.table(chr_arm_cna_path,sep="\t",header=T)
  cohort_chr_arm_cna =rbind(cohort_chr_arm_cna,chr_arm_cna)
  
#  fga_per_chrom =read.table(fga_per_chrom_path,sep="\t",header=T)
#  cohort_fga_per_chrom =rbind(cohort_fga_per_chrom,fga_per_chrom)
}
#remove _acen gvar stalk
cohort_chr_arm_cna = cohort_chr_arm_cna %>% filter(!grepl("_",chr_arm) & !grepl("X|Y",chr_arm))
#add seqnames and width
cohort_chr_arm_cna = cohort_chr_arm_cna %>% left_join(chr_arms_df[,c("chr_arm","seqnames","width")])

cohort_chr_arm_cna = cohort_chr_arm_cna  %>% call_cn_stability() %>% call_cn_alteration() %>% call_cn_copies()
#cohort_chr_arm_cna$patient_id %>% unique() %>% length()

#write.table(cohort_chr_arm_cna, cohort_chr_arm_cna_path,quote=F,sep="\t",row.names = F,col.names = T)
}


## annotate genes_cna_expression with chromosome arm info ----
#prevent auto rerun
if(FALSE){
#genes_cna_expression = read.table(genes_cna_expression_path,sep="\t",header = T)
genes_cna_expression = genes_cna_expression %>% mutate(chr_arm = ifelse(grepl("p",cytoband),paste0(seqnames,"p"),paste0(seqnames,"q"))) 
  
  chr_arm_cols = c("call","cr_stable","alteration","tumor_af")
  cohort_chr_arm_cna = cohort_chr_arm_cna %>% dplyr::rename_with(.cols=chr_arm_cols,function(x){paste0("chr_arm_",x)})
  
  genes_cna_expression = genes_cna_expression %>% left_join(cohort_chr_arm_cna[,c("patient_id","chr_arm",paste0("chr_arm_",chr_arm_cols))], by=c("patient_id","chr_arm"))
  
  ## add flag ... then fisher to compare 
  ## chr arm call == call and stable then add flag & not sv 
  genes_cna_expression = genes_cna_expression %>% mutate(helper_chr_arm_alteration_or_sv = 
                                                           ifelse( (flag_sv_bp==T | flag_sv_overlap==T), "sv",
                                                                   ifelse(call==chr_arm_call & chr_arm_cr_stable==T & call %in% c("loss","gain","loh"), "chr_arm",
                                                                          ifelse(call=="0" & flag_sv_bp==F & flag_sv_overlap==F, "none", "other" ))))              
#  write.table(genes_cna_expression,genes_cna_expression_path,sep="\t",col.names = T,row.names = F)
  }
#last step of genes cna expression annotation?

genes_cna_expression = read.table(genes_cna_expression_path,sep="\t",header = T)
  
  
# Gene centric df ----
if(length(Sys.glob(patient_gene_overview_path))==1) {
  gene_level_overview = read.table(patient_gene_overview_path,sep="\t",header=T)
} else {
gene_id_map = gene_properties_df[,c("gene_name","gene_id")] %>% unique()
gene_id_map = gene_id_map %>% filter(!gene_name %in% gene_id_map[duplicated(gene_id_map$gene_name),c("gene_name")])

snvs = cohort_snv_somatic %>% left_join(gene_id_map) %>% filter(gene_name!=""&!is.na(gene_name)) %>% 
  mutate(modality="snv") %>%
  dplyr::select(patient_id,gene_name,gene_id,modality)

cnas = cohort_gene_cna %>%  filter(call !="0") %>% 
  mutate(modality="cna")  %>%
  dplyr::select(patient_id,gene_name,gene_id,modality,call) 

svs_genes = sv_gene_long %>% filter(sv_position != "spanning") %>% mutate(modality=paste0("sv_bp")) %>% as.data.frame() %>% dplyr::select(patient_id,gene_name,gene_id,modality) 
#svs_genes_1mb = svs_1mb_gene %>% mutate(modality=paste0("sv_bp_1mb")) %>% dplyr::select(patient_id,gene_name,gene_id,modality) 

patient_gene_overview = rbind(snvs,svs_genes)#,svs_genes_1mb)

patient_gene_overview = patient_gene_overview %>% group_by(patient_id,gene_name) %>% summarize(alteration=toString(unique(sort(modality))))
patient_gene_overview = patient_gene_overview %>% full_join(cnas[,c("patient_id","gene_name","call","gene_id")])

patient_gene_overview = patient_gene_overview %>% left_join(gene_onco_tsg_consensus[,c("gene_name","onco","tsg","onco_or_tsg")]) %>% 
  left_join(gene_properties_df[,c(gene_cols,gene_properties_df_flags)]) %>% left_join(cohort[,cohort_anno_display_cols])  

patient_gene_overview = patient_gene_overview %>% annotate_gene_group(cancer_genes = cancer_genes,wilms_genes_of_interest = wilms_genes_of_interest)

write.table(patient_gene_overview,patient_gene_overview_path,sep="\t",quote=F,row.names = F,col.names = T) 
}

if(length(Sys.glob(patient_gene_overview_detailed_path))==1) {
  patient_gene_overview_detailed = read.table(patient_gene_overview_detailed_path,sep="\t",header=T)
} else {
## Detailed version will be wide dataframe for manual processing
  #checked, protein consequence also good for reverse strand genes
snvs = cohort_snv_somatic %>% filter(gene_name!=""&!is.na(gene_name)) %>% 
  mutate(protein_consequence=ifelse(Consequence=="missense_variant",
                                    paste0("p.",str_replace(Amino_acids,"/",Protein_position)),
                                    Consequence)) %>%
  group_by(patient_id,gene_name) %>% 
  summarize(snv= toString(paste0(tumor_AF," taf (",protein_consequence," ",polyphen_label,"/",sift_label,")")),
            snv_cnt=n())

genes_sv_bp = sv_gene_long %>% filter(sv_position != "spanning") %>%
  group_by(patient_id,gene_name) %>%
  summarize(sv_bp = toString(paste0(svtype," ",round(tumor_af_mean,2)," taf ",sv_position," ",
  ifelse(is.na(svlen_mean),"",
         ifelse(svlen_mean>1e6,paste0(round(svlen_mean/1e6,2)," Mbp"),
                paste0(round(svlen_mean/1e3,2)," kbp"))))),
            sv_bp_cnt = n())

#NB: 
# sv_gene_long == filtered_cohort_sv_gene_overlaps but summarized
# svs_1mb_gene == flank_gene_overlaps anti joined with  filtered_cohort_sv_gene_overlaps to remove intersects so only 1mb nearby

genes_sv_spanning = sv_gene_long %>% filter(sv_position == "spanning") %>%
  group_by(patient_id,gene_name) %>%
  summarize(sv_spanning = toString(paste0(svtype," ",round(tumor_af_mean,2)," taf ",
                                    ifelse(is.na(svlen_mean),"",
                                           ifelse(svlen_mean>1e6,paste0(round(svlen_mean/1e6,2)," Mbp"),
                                                  paste0(round(svlen_mean/1e3,2)," kbp"))))),
            sv_spanning_cnt = n())

genes_sv_1mb = summary_svs_1mb_gene %>%
  group_by(patient_id,gene_name) %>%
  summarize(sv_1mb = toString(paste0(svtype," ",round(tumor_af_mean,2)," taf dist.",round(distance_kbp_mean,2)," kbp")),
            sv_1mb_cnt = n())


patient_gene_overview_detailed = cohort_gene_cna %>% 
  dplyr::select(patient_id,gene_name,ensembl_id,call) %>% 
  left_join(genes_cna_expression %>% dplyr::select(patient_id,gene_name,ensembl_id,fpkm,vst_zscore,neutral_vst_zscore,contains("wilcox")))  %>% 
  left_join(snvs) %>%
  left_join(genes_sv_bp) %>% 
  left_join(genes_sv_spanning) %>% 
  left_join(genes_sv_1mb)

patient_gene_overview_detailed = patient_gene_overview_detailed %>%
  left_join(gene_onco_tsg_consensus[,c("gene_name","onco","tsg","onco_or_tsg")]) %>% 
  left_join(gene_properties_df[,c(gene_cols,gene_properties_df_flags)]) %>% left_join(cohort[,c("sex",cohort_anno_display_cols)])  



#relevant if  expression altered &&  sv spanning, 1mb sv, loss or gain 
#always affected if snv, sv bp, > or expression altered and then look at these cols for the cause

#need the or cohort level vst because not always CN altered 

patient_gene_overview_detailed = patient_gene_overview_detailed %>% 
  mutate(flag_expression_altered = (abs(vst_zscore)>vst_zscore_threshold | abs(neutral_vst_zscore)>vst_zscore_threshold),
         #| flag_any_wilcox),
         flag_expression_altered = ifelse(is.na(flag_expression_altered),F,flag_expression_altered))


patient_gene_overview_detailed = patient_gene_overview_detailed %>%
  dplyr::mutate(across(.cols=contains("_cnt"),.fns = function(x){ifelse(!is.na(x),x,0)})) 

## sv affected iff sv bp || (expression altered &&  sv spanning | 1mb sv)
## snv affected if snv
## cn affected if call!=0 && expression altered
## although high AF can also indicate both alleles hit > then often LOH observed 

## maybe make it multi-modal hit? or just use it manually!
patient_gene_overview_detailed = patient_gene_overview_detailed  %>%
  mutate(flag_affected_by_sv =  (sv_bp_cnt>0 | (flag_expression_altered & (sv_spanning_cnt>0 | sv_1mb_cnt>0))),
         flag_affected_by_cna = (call %in% c("gain","loss") & flag_expression_altered),
         flag_affected_by_snv =  (snv_cnt>0),
         flag_affected = (flag_affected_by_cna|flag_affected_by_sv|flag_affected_by_snv))

#potential misses because of threshold?
if(FALSE) {
patient_gene_overview_detailed %>% filter(flag_gene_of_interest & !flag_expression_altered & !flag_affected & abs(neutral_vst_zscore)>1.5) %>%
  filter(sv_bp_cnt>0 | sv_spanning_cnt>0 | sv_1mb_cnt>0 | call %in% c("gain","loss")) %>% 
  filter(flag_wilms_gene_of_interest) %>% View()
  select(patient_id,gene_name_display) 
}

patient_gene_overview_detailed$alteration=""
patient_gene_overview_detailed$alteration_simple=""

patient_gene_overview_detailed = patient_gene_overview_detailed %>% 
  mutate(alteration = ifelse(snv_cnt>0,paste0(alteration, " snv"),alteration),
         alteration = ifelse(sv_bp_cnt>0,paste0(alteration, " sv_bp"),alteration),
         alteration = ifelse((flag_expression_altered & (sv_spanning_cnt>0 | sv_1mb_cnt>0)),paste0(alteration, " sv_indirect"),alteration),
         alteration = ifelse(flag_affected_by_cna,paste0(alteration, " ",call),alteration)) %>%
  
  mutate(alteration_simple = ifelse(snv_cnt>0,paste0(alteration_simple, " snv"),alteration_simple),
         alteration_simple = ifelse(flag_affected_by_sv,paste0(alteration_simple, " sv"),alteration_simple),
         alteration_simple = ifelse(flag_affected_by_cna,paste0(alteration_simple, " cna"),alteration_simple)) %>% 
  mutate(expression_change = ifelse(flag_expression_altered,ifelse(neutral_vst_zscore>0,"up","down"),NA))

#patient_gene_overview_detailed %>% filter(gene_name=="WT1") %>% arrange(flag_affected) %>% select(patient_id,gene_name,call,neutral_vst_zscore,alteration,snv,sv_bp,sv_1mb,sv_spanning,expression_change) %>% View()

patient_gene_overview_detailed[patient_gene_overview_detailed$alteration=="",c("alteration")]=NA
patient_gene_overview_detailed[patient_gene_overview_detailed$alteration_simple=="",c("alteration_simple")]=NA



#patient_gene_overview_detailed %>% filter(gene_name=="CTNNB1" & flag_affected)
write.table(patient_gene_overview_detailed,patient_gene_overview_detailed_path,sep="\t",quote=F,row.names = F,col.names = T) 

#save a cancer genes only version for faster loading

patient_gene_overview_detailed_cancer = patient_gene_overview_detailed %>% filter(flag_gene_of_interest)

write.table(patient_gene_overview_detailed_cancer,patient_gene_overview_detailed_cancer_path,sep="\t",quote=F,row.names = F,col.names = T) 

}


## Mutation burden ----


### aneuploidy score ----

make_chromosome_summary_aneuploidy_score = function(cohort_chr_arm_cna) {

  #39 is good from some chr we only have 1 arm 
  #if both chromosome arms are same direction (gain-gain, loss-loss) then it is a whole chr gain/loss and single event. 
  
  chromosome_stable_states = c("gain_full","loss_full","cn_loh_full","iso_q","iso_p")
  chr_arm_stable_states = c("gain_stable","loss_stable","cn_loh")
  
  acrocentric_chr=c("chr13","chr14","chr15","chr21","chr22")
  
  
  
  chromosome_stable_states = c("gain_full","loss_full","cn_loh_full","iso_q","iso_p")
  chr_arm_stable_states = c("gain_stable","loss_stable","cn_loh")
  acrocentric_chr=c("chr13","chr14","chr15","chr21","chr22")
  
  
  cohort_chr_arm_cna_summary = cohort_chr_arm_cna %>% dplyr::select(seqnames,chr_arm,patient_id,alteration) %>% unique() %>%
    group_by(patient_id) %>% 
    summarize(chr_arm_cnt=n(),
              chr_arm_gains=sum(alteration=="gain_stable"),
              chr_arm_loss=sum(alteration=="loss_stable"),
              chr_arm_loh=sum(alteration=="cn_loh"))
  
  cohort_chr_arm_cna = cohort_chr_arm_cna %>%  mutate(chr_arm_orientation = ifelse(grepl("p",chr_arm),"p","q"))
  
  cohort_chromosomes = cohort_chr_arm_cna %>% dplyr::select(patient_id,seqnames,alteration,chr_arm_orientation) %>% unique() %>% pivot_wider(names_from = chr_arm_orientation,values_from=alteration)
  
  
  cohort_chromosomes = cohort_chromosomes %>% mutate(chromosome_level_event= 
                                                       ifelse((is.na(p) & seqnames %in% acrocentric_chr) | p==q, "same",
                                                              ifelse(!is.na(p)&p=="loss_stable"&q=="gain_stable","iso_q",
                                                                     ifelse(!is.na(p)&p=="gain_stable"&q=="loss_stable","iso_p",
                                                                            "unknown")))) %>% 
    mutate(alteration =
             ifelse(chromosome_level_event=="same" & q=="gain_stable","gain_full",
                    ifelse(chromosome_level_event=="same" & q=="loss_stable","loss_full",
                           ifelse(chromosome_level_event=="same" & q=="cn_loh","cn_loh_full",
                                  ifelse(chromosome_level_event=="same" & q=="0","cn_stable_full",chromosome_level_event)))))
  
  
  cohort_chromosomes_summary = 
    cohort_chromosomes %>% dplyr::select(seqnames,patient_id,alteration) %>% unique() %>%
    group_by(patient_id) %>% 
    summarize(chr_cnt=n(),
              chr_full_gain=sum(alteration=="gain_full"),
              chr_full_loss=sum(alteration=="loss_full"),
              chr_full_loh=sum(alteration=="cn_loh_full"),
              chr_full_cn_stable=sum(alteration=="cn_stable_full"),
              isochr=sum(grepl("iso",alteration)),
              chr_other=chr_cnt-(sum(grepl("full",alteration))+sum(grepl("iso",alteration))),
              chr_cnt_check=chr_other+(sum(grepl("full",alteration))+sum(grepl("iso",alteration))))
  
  # add specific chr
  
  cohort_chromosomes_wide = cohort_chromosomes %>%
    group_by(patient_id,alteration) %>% 
    summarize(chr_lst=toString(sort(unique(seqnames)))) 
  
  cohort_chromosomes_summary = cohort_chromosomes_summary %>% left_join(
    cohort_chromosomes_wide %>% pivot_wider(names_from = alteration,values_from=chr_lst) %>% dplyr::select(patient_id,gain_full,loss_full,iso_q)  %>% dplyr::rename(chr_gain_lst =gain_full,chr_loss_lst=loss_full,chr_iso_q_lst=iso_q)
  ) 
  
  cohort_chromosomes_summary 
  
  
  #For event counts: full gains, losses and isochr are single events. exclude these from chr arm level events
  
  
  cohort_chr_arm_cna_without_full = cohort_chr_arm_cna %>% dplyr::select(seqnames,chr_arm,patient_id,alteration) %>% unique() %>%
    anti_join( filter(cohort_chromosomes,alteration %in% chromosome_stable_states), by=c("patient_id","seqnames")) 
  
  cohort_chr_arm_cna_summary = cohort_chr_arm_cna_without_full %>% 
    group_by(patient_id) %>% 
    summarize(chr_arm_cnt_check=n(),
              chr_arm_gain=sum(alteration=="gain_stable"),
              chr_arm_loss=sum(alteration=="loss_stable"),
              chr_arm_loh=sum(alteration=="cn_loh"))
  
  cohort_chr_arm_cna_without_full_wide = 
    cohort_chr_arm_cna_without_full %>% 
    group_by(patient_id,alteration) %>% 
    summarize(chr_arm_lst=toString(sort(unique(chr_arm)))) 
  
  
  cohort_chr_arm_cna_summary = cohort_chr_arm_cna_summary %>% left_join(
    cohort_chr_arm_cna_without_full_wide %>% filter(alteration %in% chr_arm_stable_states) %>%  pivot_wider(names_from = alteration,values_from=chr_arm_lst) %>% dplyr::rename(chr_arm_gain_lst =gain_stable, chr_arm_loss_lst=loss_stable,chr_arm_loh_lst=cn_loh)
  )
  
  #combine
  cohort_chromosomes_summary = cohort_chromosomes_summary %>% left_join(cohort_chr_arm_cna_summary)
  
  cohort_chromosomes_summary = cohort_chromosomes_summary %>% 
    mutate(aneuploidy_score = chr_full_gain+chr_full_loss+isochr+chr_arm_gain+chr_arm_loss,
         aneuploidy_score_loh = chr_full_gain+chr_full_loss+chr_full_loh+isochr+chr_arm_gain+chr_arm_loss+chr_arm_loh)

  return(cohort_chromosomes_summary)
}

cohort_chromosomes_summary = make_chromosome_summary_aneuploidy_score(cohort_chr_arm_cna)

#write.table(cohort_chromosomes_summary,"~/results/wilms_v3/chromosome_burden.aneuploidy_score.tsv",sep="\t",row.names = F,col.names = T)

cohort_chromosomes_summary = cohort_chromosomes_summary %>% left_join(patient_id_to_labels) %>% arrange(patient_label)

cohort_chromosomes_summary_suptable = cohort_chromosomes_summary  %>% select(patient_label,aneuploidy_score_loh,contains("lst"))

#checked and is same as previously 

### combine to overview dataframe long 

cohort_triple_cols=c("patient_id","tumor_id","normal_id","rna_id")

cohort_snv_burden_somatic = read.table("~/results/wilms_v3/cohort_snv_burden_somatic.tsv",sep="\t",header = T)
ctx_burden = read.table("~/results/wilms_v3/ctx_burden.tsv",sep="\t",header = T)
export_sv_event_long = read.table("~/results/wilms_v3/sv_burden_merged_event.tsv",sep="\t",header = T)
cohort_chromosomes_summary = read.table("~/results/wilms_v3/chromosome_burden.aneuploidy_score.tsv",sep="\t",header = T)


cohort_overview = cohort_chromosomes_summary %>% select(patient_id,aneuploidy_score,aneuploidy_score_loh,contains("lst"))

cohort_overview = cohort_overview  %>% left_join(cohort_snv_burden_somatic,by = "patient_id")
cohort_overview = cohort_overview  %>% left_join(ctx_burden,by = "patient_id")

export_sv_event_long = export_sv_event_long %>% filter(cnt_type=="merged_event" & tool=="all_cnt" & svtype=="all") %>%
  select(patient_id,variable_name,cnt) %>% pivot_wider(names_from="variable_name",values_from="cnt")

cohort_overview = cohort_overview  %>% left_join(export_sv_event_long)
cohort_overview[is.na(cohort_overview)]=0

#write.table(cohort_overview,cohort_mutation_overview_path,sep = "\t",quote=F,row.names = F,col.names = T)

display_cnt_type_cols = c("tmb","snv_indel_somatic_nonsyn_cnt","somatic_sv_taf_cnt","aneuploidy_score_loh","ctx_chrom_cnt","links_short")

display_cohort_overview_wide =  cohort_overview  %>% left_join(patient_id_to_labels)  %>%
  dplyr::select(patient_label,display_cnt_type_cols,contains("lst")) %>% arrange(patient_label)
display_cohort_overview_wide$tmb = round(display_cohort_overview_wide$tmb,2)
write.table(display_cohort_overview_wide,mutation_burden_sup_table_path,sep = "\t",quote=F,row.names = F,col.names = T)


# Wnt pathway score ----
library(msigdbr)
msigdb_c2  = msigdbr::msigdbr(species="Homo sapiens",category="C2")
wnt_pathway = msigdb_c2 %>% filter(gs_id=="M39669")
wnt_pathway_genes = wnt_pathway %>% dplyr::rename(gene_name=gene_symbol, ensembl_id = ensembl_gene)
wnt_pathway_genes = wnt_pathway_genes %>% left_join(gene_properties_df[,gene_cols])
wnt_pathway_genes$gene_name_display2 = paste0(wnt_pathway_genes$gene_name,"_",wnt_pathway_genes$cytoband)

genes_vst_wnt_all = subset_vst(genes_vst,wnt_pathway_genes,c("ensembl_id","gene_name_display2"))

wnt_all_mean = colMeans(genes_vst_wnt_all)
wnt_all_mean =wnt_all_mean %>% as.data.frame()
names(wnt_all_mean)=c("colmeans_wnt_all")
wnt_all_mean$patient_id=rownames(wnt_all_mean)
wnt_all_mean = wnt_all_mean %>% left_join(patient_id_to_labels)

wnt_all_mean = wnt_all_mean %>% mutate(wnt_score = rescale(wnt_all_mean$colmeans_wnt_all, to=c(0,1)))

#write.table(wnt_all_mean,"~/results/wilms_v3/wnt_score.tsv",sep="\t",col.names = T,row.names = F,quote = F)




# Table exports for manuscript related to expression and mutations  ----

  ## chromosomes export  for manuscript ----
  cohort_chr_arm_cna = read.table(cohort_chr_arm_cna_path,sep="\t",header = T)
  
  table_cohort_chr_arm_cna = cohort_chr_arm_cna %>% left_join(patient_id_to_labels)
  
  
  table_cohort_chr_arm_cols = c("patient_label","chr_arm","cr_l2fc_50","maf_50","alteration","nr_of_copies","tumor_af")
  
  table_cohort_chr_arm_cna = table_cohort_chr_arm_cna[,table_cohort_chr_arm_cols]
  table_cohort_chr_arm_cna$chr_arm = factor(table_cohort_chr_arm_cna$chr_arm,
                                            levels=chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")])
  table_cohort_chr_arm_cna = table_cohort_chr_arm_cna %>% filter(alteration!="0" & alteration!="unknown") %>% arrange(patient_label, chr_arm)
  
  table_cohort_chr_arm_cna = table_cohort_chr_arm_cna %>% dplyr::rename(cr_l2fc=cr_l2fc_50,maf=maf_50,estimated_taf=tumor_af)
  #table_cohort_chr_arm_cna %>% View()
  write.table(table_cohort_chr_arm_cna, chr_arm_sup_table_path,quote=F,sep="\t",row.names = F,col.names = T)
  
  
  wilcox_significant_all = read.table(wilcox_overall_path,sep="\t",header = T)
  
  #export for paper, gain and loss only 
  wilcox_significant_suptable = wilcox_significant_all[!grepl("sv",wilcox_significant_all$analysis),
                                                       c("analysis","ensembl_id","gene_name","cytoband","flag_gene_of_interest","l2fc_mean","l2fc_median","pvalue","BH"
                                                       )] %>% unique() 
  
  wilcox_significant_suptable = wilcox_significant_suptable %>% dplyr::rename(FDR=BH)
  write.table(wilcox_significant_suptable,wilcox_significant_sup_table_path,sep="\t",col.names = T,row.names = F,quote = F)
  
  
  
  expr_cluster_dir="~/results/wilms_v3/expression_cluster/"
  top_n=1000
  ontologies = c("BP","CC","MF")
  
  ontology="BP"
  #for(ontology in ontologies) {
  compare_cluster_df_path = paste0(expr_cluster_dir,"metagenes_top_",top_n,".",ontology,".tsv")
  compare_cluster_results = read.table(compare_cluster_df_path,sep="\t",header = T)
  compare_cluster_results = compare_cluster_results %>% dplyr::rename(expression_cluster=Cluster)
  compare_cluster_results$expression_cluster= str_replace(compare_cluster_results$expression_cluster,"factor_","EX")
  
  compare_cluster_results= compare_cluster_results %>% arrange(expression_cluster, qvalue)
  
  write.table(compare_cluster_results, metagenes_bp_sup_table_path,quote=F,sep="\t",row.names = F,col.names = T)
  
  #}
  
  
  cn_sv_expression_changes_path = paste0(expr_cluster_dir,"cn_sv_mediated_expression_changes.",cohort_label,".tsv")
  metagene_intersect = read.table(cn_sv_expression_changes_path,sep="\t",header=T)
  metagene_intersect = metagene_intersect %>% filter(patient_cnt>1)
  metagene_intersect = metagene_intersect %>% left_join(patient_id_to_labels)
  
  metagene_intersect$expression_cluster= str_replace(metagene_intersect$expression_cluster,"factor_","EX")
  metagene_intersect=metagene_intersect[,c("expression_cluster","ensembl_id","gene_name","cytoband","patient_label","call","flag_sv_bp","neutral_vst_zscore")] %>% unique() 
  metagene_intersect$gene_name %>% unique() %>% length()
  metagene_intersect = metagene_intersect %>% arrange(expression_cluster,gene_name,patient_label)
  
  metagene_intersect = metagene_intersect %>% dplyr::rename(cna=call)
  
  
  write.table(metagene_intersect, recurrent_cn_sv_overexpressed_sup_table_path,quote=F,sep="\t",row.names = F,col.names = T)
  
  
  chr1q_overexpressed_genes_1q_tumors_cancer_path = "~/results/wilms_v3/chr1q_overexpressed_genes_1q_tumors_cancer.tsv"
  chr1q_overexpressed_genes_1q_tumors_cancer = read.table(chr1q_overexpressed_genes_1q_tumors_cancer_path,sep="\t",header=T)
  chr1q_overexpressed_genes_1q_tumors_cancer = chr1q_overexpressed_genes_1q_tumors_cancer %>% left_join(patient_id_to_labels) %>% arrange(patient_label)
  chr1q_overexpressed_genes_1q_tumors_cancer = chr1q_overexpressed_genes_1q_tumors_cancer %>% dplyr::select(patient_label,gene_name,ensembl_id,fpkm,neutral_vst_zscore,vst_zscore)
  
  write.table(chr1q_overexpressed_genes_1q_tumors_cancer,chr1q_overexpressed_genes_1q_tumors_cancer_sup_table_path,sep="\t",quote=F,row.names = F,col.names = T) 
  
  
  patient_gene_overview_detailed_cancer_sup_table = read.table(patient_gene_overview_detailed_cancer_path,sep="\t",header=T)
  
  cols_patient_gene_overview_detailed_cancer = names(patient_gene_overview_detailed_cancer)
  cols_patient_gene_overview_detailed_cancer = cols_patient_gene_overview_detailed_cancer[!grepl("wilcox",cols_patient_gene_overview_detailed_cancer)]
  cols_patient_gene_overview_detailed_cancer = cols_patient_gene_overview_detailed_cancer[!cols_patient_gene_overview_detailed_cancer %in% c("sex",cohort_anno_display_cols)]
  cols_patient_gene_overview_detailed_cancer = cols_patient_gene_overview_detailed_cancer[!cols_patient_gene_overview_detailed_cancer %in% c("onco","tsg","onco_or_tsg")]
  cols_patient_gene_overview_detailed_cancer = cols_patient_gene_overview_detailed_cancer[!cols_patient_gene_overview_detailed_cancer %in% c("gene_type","gene_coord","seqnames","gene_name_display")]
  
  patient_gene_overview_detailed_cancer_sup_table = patient_gene_overview_detailed_cancer_sup_table %>% left_join(patient_id_to_labels)
  
  patient_gene_overview_detailed_cancer_sup_table = patient_gene_overview_detailed_cancer_sup_table %>% filter(flag_affected)
  patient_gene_overview_detailed_cancer_sup_table = patient_gene_overview_detailed_cancer_sup_table[,c("patient_label",cols_patient_gene_overview_detailed_cancer)] %>% arrange(patient_label)
  
  patient_gene_overview_detailed_cancer_sup_table = patient_gene_overview_detailed_cancer_sup_table %>% dplyr::rename(cna=call)
  
  write.table(patient_gene_overview_detailed_cancer_sup_table,patient_gene_overview_detailed_cancer_sup_table_path,sep="\t",quote=F,row.names = F,col.names = T) 

  #wnt pathway enrichment => not used
  wnt_gene_wilcoxon = read.table("~/results/wilms_v3/expression_wilcox.wnt_pathway_genes.pretreated_bm_outlier.tsv",sep="\t",header=T)
  
  export_wnt_gene_ex1 = wnt_gene_wilcoxon %>% filter(BH<fdr_wilcox_threshold & abs(l2fc_median) > l2fc_wilcox_threshold) %>% arrange(-l2fc_median)
  
  write.table(export_wnt_gene_ex1,wnt_pathway_ex1_sup_table_path,sep="\t",quote=F,row.names = F,col.names = T) 
  
  
  #from mutation burden Rmd
  kruskal_tests_df = read.table("~/results/wilms_v3/mutation_burden_overview.kruskal_tests.tsv",sep="\t",header=T)
  
  summary_stats_df = read.table("~/results/wilms_v3/mutation_burden_overview.summary_stats.tsv",sep="\t",header=T)
  
  stats_table = kruskal_tests_df %>% left_join(summary_stats_df,by = c("mutation_type"="cnt_type","variable"))
  
  stats_table = stats_table %>% mutate(value=str_replace(value,"factor_","EX"))
  stats_table = stats_table %>% mutate(variable=str_replace(variable,"tumor_stage_tests","tumor_stage"),
                                       variable=str_replace(variable,"cn_similarity_cluster_display","cn_similarity_cluster"),
                                       variable=str_replace(variable,"histology_consensus","histology"))
  
  mutation_type_reporting_lst = c("snv_indel_somatic_nonsyn_cnt","somatic_sv_taf_cnt","aneuploidy_score_loh","fga")
  variable_reporting_lst = c("histology","risk_group","tumor_stage")
  #"expression_cluster","cn_similarity_cluster",
  
  stats_table = stats_table %>% filter(mutation_type %in% mutation_type_reporting_lst & variable %in% variable_reporting_lst)
  
  stats_table = stats_table %>% dplyr::rename(kruskal_pvalue=k.pval)
  
  write.table(stats_table,mutation_burden_statistics_sup_table_path,quote=F,sep="\t",row.names = F,col.names = T)
  
  
  
  ## wnt supplementary table
  wnt_supplementary_table = patient_gene_overview_detailed %>% filter(flag_affected) %>% filter(gene_name %in% c("CTNNB1","WT1","AMER1"))

  wnt_supplementary_table = wnt_supplementary_table  %>% group_by(patient_id,gene_name) %>% 
    summarize(alteration="",
              alteration = ifelse(!is.na(snv),paste0(alteration, " snv:",snv),alteration),
              alteration = ifelse(!is.na(sv_bp),paste0(alteration, "; sv_bp:",sv_bp),alteration),
              alteration = ifelse(!is.na(sv_spanning),paste0(alteration, "; sv_spanning:",sv_spanning),alteration),
              alteration = ifelse(!is.na(sv_1mb),paste0(alteration, "; sv_1mb:",sv_1mb),alteration),
              alteration = ifelse(call!="0",paste0(alteration,"; cna:",call),alteration),
              alteration = ifelse(flag_expression_altered,paste0(alteration, " expression:",expression_change,"(",round(vst_zscore,2)," zscore, ",round(neutral_vst_zscore,2)," nzscore)"),alteration),
              alteration=trimws(alteration)) %>% 
    pivot_wider(names_from = gene_name,values_from = alteration) %>%
    left_join(cohort[,c("patient_id","patient_label","expression_cluster","c1a_pattern")]) %>% select(-patient_id) %>%
    mutate(expression_cluster=str_replace(expression_cluster,"factor_","EX")) 
  
  
  #part manual, so do not override
  # write.table(wnt_supplementary_table,"~/results/wilms_v3/wnt_supplementary_table.wilms_v3.tsv",quote=F,sep="\t",row.names = F,col.names = T)
  
  #double check:
  # M991AAA	p.S45del
  # M035AAD	p.H36P
  # M407AAA	p.S45P
  # M403AAB	p.T41A
  # M472AAC	p.S45F
  
  
  merged_svs = read.table(merged_svs_path,sep="\t",header=T)
  merged_svs = merged_svs %>% left_join(patient_id_to_labels)  %>% filter(flag_in_filtered_svs) %>% 
    filter(chrom %in% autosomes & (partner_chrom == "" |  partner_chrom %in% autosomes))
  
  sv_export_cols=c("patient_label","svtype","sv_merged_coordinate","tumor_af_mean","svlen_mean",
                   "svlen_mean_str","chrom","partner_chrom","chr_arm","partner_chr_arm",
                   "cytoband","partner_cytoband","coordinates","tools","c1a_pattern")
  

  merged_svs$chr_arm = factor(merged_svs$chr_arm,levels=chr_arms_df$chr_arm)
  merged_svs$partner_chr_arm = factor(merged_svs$partner_chr_arm,levels=chr_arms_df$chr_arm)
  
  merged_svs = merged_svs  %>% select(all_of(sv_export_cols)) %>% arrange(chr_arm) 

  merged_svs[merged_svs$c1a_pattern=="1q_gain_only",c("c1a_pattern")]='1q_gain'
  
  
  write.table(merged_svs,somatic_svs_sup_table_path,quote=F,sep="\t",row.names = F,col.names = T)


