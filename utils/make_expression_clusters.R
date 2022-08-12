library(tidyverse, quietly=TRUE)
library("matrixStats")
library(NMF)
library(DESeq2)
library(stringi)
library(GenomicRanges)


remove_version_from_id = function(identifier) {
  return(sub("\\..*","",identifier))
} 

resources_dir="/Users/ianthevanbelzen/Documents/resources/"
wdir="/Users/ianthevanbelzen/PycharmProjects/gene_expression/"
metadata_dir =  "/Users/ianthevanbelzen/data/metadata/"
reports_dir_base= "/Users/ianthevanbelzen/PycharmProjects/gene_expression/wilms/"

source("/Users/ianthevanbelzen/PycharmProjects/structuralvariation/sv_functional_analysis/functions.expression.R")
source("/Users/ianthevanbelzen/PycharmProjects/structuralvariation/sv_functional_analysis/functions.svs.R")
source("/Users/ianthevanbelzen/PycharmProjects/structuralvariation/sv_functional_analysis/functions.general.R")

source(paste0(wdir,"functions.deg_analysis.R"))

source("/Users/ianthevanbelzen/PycharmProjects/structuralvariation/sv_functional_analysis/wilms.conf")

nmf_reports_dir="/Users/ianthevanbelzen/PycharmProjects/gene_expression/wilms/nmf/"

nmf_analysis_label="nmf_exp_cluster_brunet.vst."
gene_set_label="top10000.noX"

cohort_label="pretreated_bm_outlier"
results_name="pretreated_bm_outlier_top10k_noX_blinded"


fullset_genes_vst_file = "genes.vst.tsv"
genes_vst_blinded_path = paste0(nmf_reports_dir,cohort_label,".blinded.",fullset_genes_vst_file)
genes_vst_notblinded_path=paste0(nmf_reports_dir,cohort_label,".",fullset_genes_vst_file)

nmf_reports_dir_base="/Users/ianthevanbelzen/PycharmProjects/gene_expression/wilms/nmf/"


## Read cohort ----
cohort_sup_table_path="~/data/metadata/cohort.wilms_v3.supplementary_table_1.tsv"

cohort = read.table(cohort_sup_table_path,sep = "\t", header=T) 
cohort=cohort %>% filter(flag_use_expression)
cohort$patient_id=cohort$patient_label
cohort %>% nrow()

annotation = cohort[,c("patient_id","rna_id","sex","age","cn_similarity_cluster")]
annotation$sex=factor(annotation$sex,levels= c("male","female"))
annotation$age=as.numeric(annotation$age)
annotation$cn_similarity_cluster=factor(annotation$cn_similarity_cluster)
rownames(annotation)=annotation$patient_id


### Resources ----

map_template_vars=c('${resources_dir}'=resources_dir)#,#'${data_dir}'=data_dir,
#'${merged_svs_dir}'=merged_svs_dir,'${cohort_wdir}'=cohort_wdir,
#'${cohort_identifier}'=cohort_identifier)

reference = "GRCh38_gencode_v31_CTAT_lib_Oct012019"
gtf_path_template = paste0("${resources_dir}ref_annot_",reference,".gtf.gz")
gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)

gene_properties_df = get_gene_properties_df(gtf_path)
gene_properties_df$gene_coord = gene_properties_df$to_coordinate

chromosome_bands_path_template="${resources_dir}chromosome_bands.gz"
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband

gene_properties = GRanges(gene_properties_df)
names(gene_properties) = gene_properties$gene_id


gene_to_cytoband_mapping = get_genes_to_cytobands(gene_properties,chromosome_bands)

gene_properties_df = gene_properties_df %>% left_join(gene_to_cytoband_mapping)
gene_properties_df$gene_name_display =  paste0(gene_properties_df$gene_name,"_",gene_properties_df$cytoband,"_",gene_properties_df$ensembl_id)
gene_cols=c("ensembl_id","gene_name","gene_type","gene_coord","seqnames","cytoband","gene_name_display")





## Use DESeq2 normalized cnts ----
## same but all patients > perform vst here!
### Make vst normalized readcount ----

#tried different designs, but decided on blinded 

dds_design = ~ 1

expression_df = read.table(cohort_gene_expression_path,sep="\t",header=T)
expression_overview_wide = prepare_expression_df_deseq2(expression_df,cohort$patient_id)

if(length(names(expression_overview_wide)) !=  length(cohort$patient_id)) {
  print("Mismatch expression overview and cohort")
  #quit()
}

## Make DDS object ----
dds = DESeqDataSetFromMatrix(countData = expression_overview_wide,
                             colData = annotation,
                             design = dds_design )

nrow(dds)
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep,]
nrow(dds)

dds =DESeq(dds)

## fullset genes, save normalized counts ----

#blinded
vst = varianceStabilizingTransformation(dds,blind=T)
  
fullset_genes = as.data.frame(assay(vst))
fullset_genes$ensembl_id = rownames(fullset_genes)
write.table(fullset_genes,genes_vst_blinded_path,quote=F,sep="\t",row.names = F)

fullset_genes = fullset_genes %>% left_join(gene_properties_df[!grepl("PAR_Y",gene_properties_df$gene_id),c("ensembl_id","gene_name_display")] %>% unique())

genes_vst_names = fullset_genes 
rownames(genes_vst_names) =genes_vst_names$gene_name_display
genes_vst_names = genes_vst_names %>% dplyr::select(where(is.numeric))


nmf_plot_cluster_path = paste0(nmf_reports_dir,nmf_analysis_label,gene_set_label,".",cohort_label,".blinded.pdf")
nmf_results_path= paste0(nmf_reports_dir,nmf_analysis_label,gene_set_label,".",cohort_label,".blinded.rds")


  cluster_range=2:7
  #remove y, Mt genes
  mostvar_expr_vst = get_most_variable_rows(genes_vst_names[!grepl("chrY|chrX|chrM",rownames(genes_vst_names)),],n=10000)
  
  exp_cluster_estimate_brunet <- nmf(mostvar_expr_vst, rank = cluster_range, method = "brunet", 
                                     nrun = 10, seed = 123456, .opt = "v-p")
  
  
  pdf(nmf_plot_cluster_path)
  print( plot(exp_cluster_estimate_brunet) ) 
  dev.off()
  
  saveRDS(exp_cluster_estimate_brunet,nmf_results_path)
  

  
  ## Blind = F
  
  
  # unblinded 
 
  vst = varianceStabilizingTransformation(dds,blind=FALSE)
  
  fullset_genes = as.data.frame(assay(vst))
  fullset_genes$ensembl_id = rownames(fullset_genes)
  write.table(fullset_genes,genes_vst_notblinded_path,quote=F,sep="\t",row.names = F)
  
  
  fullset_genes = fullset_genes %>% left_join(gene_properties_df[!grepl("PAR_Y",gene_properties_df$gene_id),c("ensembl_id","gene_name_display")] %>% unique())
  
  genes_vst_names = fullset_genes 
  rownames(genes_vst_names) =genes_vst_names$gene_name_display
  genes_vst_names = genes_vst_names %>% dplyr::select(where(is.numeric))
  
  
  
  nmf_plot_cluster_path = paste0(nmf_reports_dir,nmf_analysis_label,gene_set_label,".",cohort_label,".pdf")
  nmf_results_path= paste0(nmf_reports_dir,nmf_analysis_label,gene_set_label,".",cohort_label,".rds")
  
  
  cluster_range=2:7
  #remove y, Mt genes
  mostvar_expr_vst_2 = get_most_variable_rows(genes_vst_names[!grepl("chrY|chrX|chrM",rownames(genes_vst_names)),],n=10000)
  
  exp_cluster_estimate_brunet <- nmf(mostvar_expr_vst_2, rank = cluster_range, method = "brunet", 
                                     nrun = 10, seed = 123456, .opt = "v-p")
  
  pdf(nmf_plot_cluster_path)
  print( plot(exp_cluster_estimate_brunet) ) 
  dev.off()
  
  saveRDS(exp_cluster_estimate_brunet,nmf_results_path)
  
  


### Analyse results  ----


  nmf_results_path= paste0(nmf_reports_dir_base,nmf_analysis_label,gene_set_label,".",cohort_label,".blinded.rds")
  genes_vst_path= paste0(nmf_reports_dir_base,cohort_label,".blinded.",fullset_genes_vst_file)
  
  ### Read data ----
  
  nmf_reports_dir= paste0(nmf_reports_dir_base,results_name,"/")
  load_vst(genes_vst_path,as_global=T)
  
  ntop=10000
  mostvar_expr_vst = get_most_variable_rows(genes_vst_names[!grepl("chrY|chrM|chrX",rownames(genes_vst_names)),],n=ntop)
  
  
  exp_cluster_estimate_brunet= readRDS(nmf_results_path)
  print( plot(exp_cluster_estimate_brunet) ) 
  
  nmf_cluster_k=4
  nmf_fit = exp_cluster_estimate_brunet$fit
  nmf_fit = nmf_fit[names(nmf_fit) == nmf_cluster_k]
  nmf_fit = nmf_fit[[1]]
  
  
  nmf.h <- NMF::basis(nmf_fit)
  
  nmf.w <- NMF::coef(nmf_fit)
  rownames(nmf.w) = paste0("factor_",1:nrow(nmf.w))
  
  nmfw <- t(nmf.w)
  nmf_df <- as.data.frame(nmf.h)
  colnames(nmf_df) =paste0("factor_",1:length(nmf_df))
  nmf_df_plot = nmf_df - rowMeans(nmf_df)
  
  #patient vs metagene
  #placed into a cluster corresponding to the most highly expressed metagene in the sample; 
  #that is, sample j is placed in cluster i if the hij is the largest entry in column
  patients_vs_metagenes = as.data.frame(nmfw)
  patients_vs_metagenes$expression_cluster = colnames(patients_vs_metagenes)[max.col(patients_vs_metagenes, ties.method = "first")]
  patients_vs_metagenes$patient_id = rownames(patients_vs_metagenes)

  
  metagenes = as.data.frame(nmf_df)
  metagenes$expression_cluster = colnames(metagenes)[max.col(metagenes, ties.method = "first")]
  metagenes$gene_name_display = rownames(metagenes)
  
  metagenes = metagenes %>% left_join(gene_properties_df[,c("gene_name_display","ensembl_id","seqnames","cytoband")])
  
  metagenes_foldchange= data.frame()
  for(factor in unique(metagenes$expression_cluster)) {
    metagene_genes = metagenes %>% filter(expression_cluster==factor)
    
    fold_change = rowMeans(as.matrix(genes_vst_names[metagene_genes$gene_name_display,
                                                     filter(patients_vs_metagenes,expression_cluster==factor)$patient_id])) / 
      rowMeans(genes_vst_names[metagene_genes$gene_name_display,
                               filter(patients_vs_metagenes,expression_cluster!=factor)$patient_id])
    
    fold_change = as.data.frame(fold_change)
    fold_change$expression_cluster=factor
    fold_change$gene_name_display = rownames(fold_change)
    
    metagenes_foldchange = rbind(metagenes_foldchange,fold_change)
  }
  
  
  metagenes_foldchange$l2fc=log2(metagenes_foldchange$fold_change)
  
  metagenes_foldchange = metagenes_foldchange %>% mutate(direction=ifelse(l2fc>0,"up","down"))
  
  metagenes = metagenes %>% left_join(metagenes_foldchange,by=c("gene_name_display","expression_cluster"))
  
  
### Export ----
  metagenes_path = paste0(nmf_reports_dir,"metagenes.",results_name,".tsv" )
  write.table(metagenes,metagenes_path,sep="\t",row.names = F,col.names = T,quote=F)
  
  patients_vs_metagenes_path = paste0(nmf_reports_dir,"patients_vs_metagenes.",results_name,".tsv" )
  write.table(patients_vs_metagenes,patients_vs_metagenes_path,sep="\t",row.names = F,col.names = T,quote=F)

## Metagenes:GSEA and enrichments ---
  
## took from Rmd file

  
  #Functions
  
  get_top_x = function(metagenes,ex_clusters,top_n=50) {
    metagenes_top_x=data.frame()
    for(ex_cluster in ex_clusters) {
      top_metagenes = metagenes %>% filter(expression_cluster==ex_cluster) %>% arrange(-abs(l2fc)) %>% head(n=top_n)
      metagenes_top_x = rbind(metagenes_top_x,top_metagenes)
    }
    return(metagenes_top_x)
  }

  
  # Read files
  
#  Read all genes vst and get 10k most var

  
  load_vst(genes_vst_path,as_global=T)
  
  genes_vst_names = fullset_genes 
  rownames(genes_vst_names) =genes_vst_names$gene_name_display
  genes_vst_names = genes_vst_names %>% dplyr::select(where(is.numeric))
  
  #remove y, Mt genes
  mostvar_expr_vst = get_most_variable_rows(genes_vst_names[!grepl("chrY|chrX|chrM",rownames(genes_vst_names)),],n=10000)

#  Read meta genes

  metagenes = read.table(metagenes_path,sep="\t",header=T)
  metagenes$expression_cluster = factor(metagenes$expression_cluster,levels=expression_cluster_order)
  metagenes = metagenes %>% left_join(gene_properties_df[,c("gene_name_display","gene_name","cytoband")])
  metagenes$gene_name_display2 = paste0(metagenes$gene_name,"_",metagenes$cytoband)
  
  metagenes_top_x = metagenes %>% get_top_x(ex_clusters = expression_cluster_order, top_n = 50)
  
  map_metagenes_top_x = metagenes_top_x %>% group_by(ensembl_id) %>% summarize(expr_cluster = toString(sort(unique(expression_cluster))))
  map_metagenes_top_x %>% filter(grepl(",",expr_cluster) ) %>% nrow() ==0 
  map_metagenes_top_x %>% group_by(expr_cluster) %>% count()
  
  metagenes$expr_cluster=metagenes$expression_cluster
  metagenes %>% group_by(expr_cluster) %>% count()
  
  cohort %>% group_by(expression_cluster) %>% count()
  
  metagenes = metagenes %>% mutate(flag_top_50 = ensembl_id %in% map_metagenes_top_x$ensembl_id)
  
 
  library(clusterProfiler)
  library(AnnotationHub) 
  hub = AnnotationHub()
  ahDb = query(hub, c("apiens", "OrgDb"))
  orgDb=ahDb[["AH92581"]]
  
  
  library(ReactomePA)
  library(org.Hs.eg.db)
  
  
  ontologies = c("BP","CC","MF")
  
  
  #as universe use most var genes that went as input in the pipeline
  most_var_gene_df = mostvar_expr_vst
  most_var_gene_df$gene_name_display = rownames(most_var_gene_df)
  most_var_gene_df = most_var_gene_df %>% left_join(gene_properties_df[,c("gene_name_display","ensembl_id")])
  
  universe_ensembl_to_entrez = bitr(most_var_gene_df[,c("ensembl_id")], fromType="ENSEMBL", toType="ENTREZID",OrgDb = orgDb) %>% dplyr::rename(ensembl_id=ENSEMBL)
  most_var_gene_df = most_var_gene_df  %>% left_join(universe_ensembl_to_entrez)
  universe_entrez = most_var_gene_df %>% filter(!is.na(ENTREZID))
  
  universe_gene_lst = most_var_gene_df$ensembl_id
  universe_gene_lst_entrez = universe_entrez$ENTREZID
  
  
  
  
  gene_cluster = metagenes %>% left_join(universe_ensembl_to_entrez) %>% filter(!is.na(ENTREZID))
  top_n=1000
  gene_cluster = gene_cluster %>% get_top_x(top_n,ex_clusters = expression_cluster_order)
  gene_cluster = gene_cluster %>% group_by(expression_cluster) %>% summarize(entrez_lst=list(ENTREZID))
  
  gene_cluster_lists = gene_cluster$entrez_lst
  names(gene_cluster_lists) = gene_cluster$expression_cluster
  

    for(ontology in ontologies) {
      compare_cluster_rds_path = paste0(plot_dir,"metagenes_top_",top_n,".",ontology,".rds")
      compare_cluster_df_path = paste0(plot_dir,"metagenes_top_",top_n,".",ontology,".tsv")
      compare_cluster_pdf_path = paste0(plot_dir,"metagenes_top_",top_n,".",ontology,".pdf")
      
      ck <- compareCluster(geneCluster = gene_cluster_lists, fun = enrichGO,OrgDb=orgDb,ont=ontology,universe=universe_gene_lst_entrez)
      ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      
      pdf(compare_cluster_pdf_path)
      p= dotplot(ck)
      print(p)
      dev.off()
      
      saveRDS(ck,compare_cluster_rds_path)
      write.table(as.data.frame(ck),compare_cluster_df_path,quote=F,sep="\t",row.names = F)
    }
    
  
 



