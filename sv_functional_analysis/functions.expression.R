#functions.expression.R

get_deseq_annotation = function(cohort) {
  annotation = cohort %>% select(patient_id,rna_id,tumor_type,primary_group,domain,primary_group_shorthand_label,domain_label)
  rownames(annotation) = annotation$patient_id
  annotation$tumor_type = factor(annotation$tumor_type)
  annotation$primary_group = factor(annotation$primary_group)
  annotation$domain = factor(annotation$domain)
  annotation$primary_group_label = paste0(annotation$primary_group,"-",annotation$primary_group_shorthand_label)
  annotation$primary_group_label =  factor(annotation$primary_group_label, levels = unique(annotation$primary_group_label[order(annotation$primary_group)]))
  annotation$domain_label = paste0(annotation$domain,"-",annotation$domain_label)
  annotation$domain_label =  factor(annotation$domain_label, levels = unique(annotation$domain_label[order(annotation$domain)]))
  annotation = annotation[!is.na(annotation$domain),]
  return(annotation)
}

read_expression_data = function(gene_expr_filepath) {
  gene_expression_df = read.table(gene_expr_filepath,sep="\t",header = F)
  colnames(gene_expression_df) = c("ensembl_id","counts","cpm","fpkm","chr","start","end","strand","length","ensembl_id_copy","gene_name","transcript_id")
  
  ##remove  PAR_Y genes before removing version 
  gene_expression_df = gene_expression_df[grepl("ENSG",gene_expression_df$ensembl_id)&!grepl("PAR_Y",gene_expression_df$ensembl_id),]
  gene_expression_df$gene_id=gene_expression_df$ensembl_id
  gene_expression_df$ensembl_id = remove_version_from_id(gene_expression_df$ensembl_id)
  
  return(gene_expression_df)
}

get_cohort_gene_expression = function(sample_lst,expression_data_dir= "/Users/ianthevanbelzen/data/rna_expression/gencode31oct2019_multioverlap_largest_overlap/") {
  
  cohort_gene_expression = data.frame(stringsAsFactors = F)
  for(rna_id in sample_lst) {
    gene_expr_filepath = paste0(expression_data_dir,rna_id,expression_gene_file_ext)
    gene_expression_df = read_expression_data(gene_expr_filepath)
    gene_expression_df$rna_id = rna_id
    
    cohort_gene_expression = rbind(cohort_gene_expression,gene_expression_df[,c("ensembl_id","gene_id","gene_name","counts","cpm","fpkm","rna_id")])
  }
  
  cohort_gene_expression$ensembl_id=factor(cohort_gene_expression$ensembl_id)
  return(cohort_gene_expression)
}

## TODO function to remove lowly expressed across samples

#RELATIVE TO GROUP, so patient compared to other patients of that group
calculate_fpkm_zscore = function(gene_expression,group_cols=c("ensembl_id")){
  if(!"fpkm_log" %in% names(gene_expression)) {
    gene_expression = gene_expression  %>% mutate(fpkm_log = log(fpkm+0.001))
  }
  #note that this is mutate so adds group cols but still keeps the patient level data 
  grouped_gene_expression = gene_expression %>%
    group_by(across(all_of(group_cols))) %>% 
    mutate(sum_counts = sum(counts),
           fpkm_mean = mean(fpkm_log,na.rm = T), 
           fpkm_sd = sd(fpkm_log,na.rm = T), 
           fpkm_zscore = (fpkm_log - fpkm_mean )/ fpkm_sd) #%>% ungroup() %>% unique()
  
  gene_level_normal_distribution = grouped_gene_expression %>% 
    group_by(across(all_of(group_cols))) %>%  
    summarize(group_size=n(),
              shapiro_pval =  ifelse((group_size>2 & sum_counts>0),shapiro.test(fpkm_log)$p.value,0),
              fpkm_log_normal_dist=(shapiro_pval>=0.05)) %>% unique()
  
  grouped_gene_expression = grouped_gene_expression %>%
    left_join(gene_level_normal_distribution, by=group_cols) %>% unique()
  
  return(grouped_gene_expression)
}


## NOTE: not so usefull.... > do wilcox instead for group comparisons
calculate_fpkm_zscore_groupwise = function(gene_expression,group_cols=c("ensembl_id")){
  if(!"fpkm_log" %in% names(gene_expression)) {
    gene_expression = gene_expression  %>% mutate(fpkm_log = log(fpkm+0.001))
  }
  #note that this is mutate so adds group cols but still keeps the patient level data 
  gene_expression_summary = gene_expression %>%
    group_by(across(all_of(group_cols))) %>% 
    summarize(fpkm_log_group = mean(fpkm_log,na.rm = T),.groups="drop")
  
  compare_gene_expression = gene_expression_summary %>%
    group_by(ensembl_id) %>% 
    mutate(
           fpkm_mean = mean(fpkm_log_group,na.rm = T), 
           fpkm_sd = sd(fpkm_log_group,na.rm = T), 
           fpkm_zscore = (fpkm_log_group - fpkm_mean )/ fpkm_sd) #%>% ungroup() %>% unique()
  
  return(compare_gene_expression)
}

calculate_vst_zscore = function(gene_expression,group_cols=c("ensembl_id"),attr_prefix=NULL){
  #implement as check for vst
  incomplete_genes_lst = gene_expression[is.na(gene_expression$vst),c("ensembl_id")]
  incomplete_genes_expression = gene_expression %>% filter(ensembl_id %in% incomplete_genes_lst) 
  
  gene_expression = gene_expression %>% filter(!ensembl_id %in% incomplete_genes_lst) 
  
  #note that this is mutate so adds group cols but still keeps the patient level data 
  grouped_gene_expression = gene_expression %>%
    group_by(across(all_of(group_cols))) %>% 
    mutate(sum_counts = sum(counts),
           vst_mean = mean(vst,na.rm = T), 
           vst_sd = sd(vst,na.rm = T), 
           vst_zscore = (vst - vst_mean )/ vst_sd) #%>% ungroup() %>% unique()
  
  gene_level_normal_distribution = grouped_gene_expression %>% 
    group_by(across(all_of(group_cols))) %>%  
    summarize(group_size=n(),
              shapiro_pval =  ifelse((group_size>2 & sum_counts>0 & !is.na(vst)),shapiro.test(vst)$p.value,0),
              vst_normal_dist=(shapiro_pval>=0.05)) %>% unique()
  
  grouped_gene_expression = grouped_gene_expression %>%
    left_join(gene_level_normal_distribution, by=group_cols) %>% unique()
  
  
  if(nrow(incomplete_genes_expression)>0){
    grouped_gene_expression = rbind_no_colmatch(grouped_gene_expression,incomplete_genes_expression)
  }
  
  if(!is.null(attr_prefix)) {
    new_cols = names(grouped_gene_expression)[!names(grouped_gene_expression) %in% names(gene_expression)]
    grouped_gene_expression = grouped_gene_expression %>% dplyr::rename_with(.cols=all_of(new_cols),.fn=function(x){paste0(attr_prefix,x)}) 
  }
  
  return(grouped_gene_expression)
}



call_expression = function(genes_cna_expression, fpkm_zscore_threshold=0.96){
  genes_cna_expression=genes_cna_expression %>% mutate(call_expr = 
                                                         ifelse(fpkm_zscore>fpkm_zscore_threshold,"gain",
                                                                ifelse(fpkm_zscore<(-fpkm_zscore_threshold),"loss",
                                                                       ifelse(abs(fpkm_zscore)<fpkm_zscore_threshold,"0",NA))))
  return(genes_cna_expression)
}

wilcox_test_groups_fpkm_cna = function(change_group,neutral_group,prefix=NULL){
  entry=c()
  
  wilcox=wilcox.test(change_group$fpkm,
                     neutral_group$fpkm)
  
  entry$pvalue = wilcox$p.value
  
  entry$change_group_size = change_group %>% nrow()
  entry$neutral_group_size = neutral_group %>% nrow()
  
  entry$change_fpkm_mean = change_group$fpkm %>% mean()
  entry$neutral_fpkm_mean = neutral_group$fpkm %>% mean()
  entry$change_fpkm_median = change_group$fpkm %>% median()
  entry$neutral_fpkm_median = neutral_group$fpkm %>% median()
  
  if("cr_l2fc_50" %in% names(entry)) {
    entry$change_cr_l2fc_50_mean = change_group$cr_l2fc_50 %>% mean()
    entry$neutral_cr_l2fc_50_mean = neutral_group$cr_l2fc_50 %>% mean()
    entry$change_cr_l2fc_50_median = change_group$cr_l2fc_50 %>% median()
    entry$neutral_cr_l2fc_50_median = neutral_group$cr_l2fc_50 %>% median()
  }
  if(!is.null(prefix)){
    entry_names=names(entry)
    names(entry) = str_replace(entry_names,"change_",prefix)
  }
  return(entry)
}

get_wilcox_expression_groups = function(cohort_gene_expression_long,target_gene_lst,gene_id_col="gene_id",
                                        attr_col="call",attr_value="gain",attr_value_contrast="0",prefix=NULL) {
  df = data.frame()
  target_gene_lst = target_gene_lst[target_gene_lst %in% cohort_gene_expression_long[,gene_id_col]]
  
  for(target_gene in target_gene_lst) {
    entry=c()
    
    target_expression = cohort_gene_expression_long[cohort_gene_expression_long[,gene_id_col]==target_gene,]
    if(nrow(target_expression)==0) { next() }
    test_group = target_expression[target_expression[,attr_col]==attr_value,]
    
    if(!is.null(attr_value_contrast)) {
      contrast_group = target_expression[target_expression[,attr_col]==attr_value_contrast,]
    } else {
      contrast_group = target_expression[target_expression[,attr_col]!=attr_value,]
    }
    
    
    if(nrow(test_group)==0|nrow(contrast_group)==0) { next() } 
    
    entry = wilcox_test_groups_fpkm_cna(test_group,contrast_group,prefix=prefix)
    entry[gene_id_col] = target_gene
    
    df = rbind(df,entry)
    
  } 
  
  if(nrow(df)==0) {
    print("ERROR: no df, check target gene list ")
    return()
  }
  
  if(is.null(prefix)){
    prefix="change_"
  }  
  
  df$change_fpkm_mean = df[,paste0(prefix,"fpkm_mean")]
  df$change_fpkm_median = df[,paste0(prefix,"fpkm_median")]
  
  df = df  %>% 
    mutate( l2fc_mean=log2(change_fpkm_mean/neutral_fpkm_mean),
            l2fc_median=log2(change_fpkm_median/neutral_fpkm_median)) 
  
  if(prefix!="change_") {
    #remove tmp
    df = df %>% dplyr::select(-change_fpkm_mean, -change_fpkm_median)
  }
  
  df = df %>% arrange(-abs(l2fc_mean))
  return(df)
}

load_vst = function(fullset_genes_path,as_global=F,sample_lst=NULL) {
  fullset_genes = read.table(fullset_genes_path,sep="\t",header=T)
  if(!is.null(sample_lst)){
    fullset_genes = fullset_genes[,names(fullset_genes) %in% c("ensembl_id","gene_name_display",sample_lst )]
  }
  
  genes_vst = fullset_genes 
  rownames(genes_vst) =genes_vst$ensembl_id
  genes_vst = genes_vst %>% dplyr::select(where(is.numeric))
  
  fullset_genes = fullset_genes %>% left_join(gene_properties_df[!grepl("PAR_Y",gene_properties_df$gene_id),c("ensembl_id","gene_name_display")] %>% unique())
  
  genes_vst_names = fullset_genes 
  rownames(genes_vst_names) =genes_vst_names$gene_name_display
  genes_vst_names = genes_vst_names %>% dplyr::select(where(is.numeric))
  
  genes_vst_plot= genes_vst - rowMeans(genes_vst)
  
  if(as_global) {
    fullset_genes <<- fullset_genes
    genes_vst <<- genes_vst
    genes_vst_names <<- genes_vst_names
    genes_vst_plot <<- genes_vst_plot
  } else {
    return(fullset_genes)
  }
}



make_gene_vst_unfiltered = function(sample_lst,cohort_gene_expression_path=NULL,expression_df=NULL,annotation=NULL) {
  library(DESeq2)
  if(is.null(cohort_gene_expression_path) & is.null(expression_df)){
    print("Need either path or df with expression data, exiting")
    return(data.frame())
  }
  if(is.null(expression_df)){
    expression_df = read.table(cohort_gene_expression_path,sep="\t",header=T)
  }
  
  #to make sure no genes are missed do not filter
  expression_overview_wide = prepare_expression_df_deseq2(expression_df,sample_lst,filter_read_counts = 0)
  
  if(length(names(expression_overview_wide)) !=  length(sample_lst)) {
    print("Mismatch expression overview and cohort")
    #quit()
  }
  
  dds_design = ~ 1
  
  if(is.null(annotation)) {
    annotation=as.data.frame(sample_lst)
    rownames(annotation) = sample_lst
  }
  ## Make DDS object
  dds = DESeqDataSetFromMatrix(countData = expression_overview_wide[sample_lst],
                               colData = annotation,
                               design = dds_design )
  
  dds =DESeq(dds)
  
  #blinded as default
  vst = varianceStabilizingTransformation(dds,blind=T)
  
  fullset_genes = as.data.frame(assay(vst))
  fullset_genes$ensembl_id = rownames(fullset_genes)
  return(fullset_genes)
}

wilcox_test_groups_vst_cna = function(change_group,neutral_group,prefix=NULL){
  entry=c()
  
  wilcox=wilcox.test(change_group$vst,
                     neutral_group$vst)
  
  entry$pvalue = wilcox$p.value
  
  entry$change_group_size = change_group %>% nrow()
  entry$neutral_group_size = neutral_group %>% nrow()
  
  entry$change_vst_mean = change_group$vst %>% mean()
  entry$neutral_vst_mean = neutral_group$vst %>% mean()
  entry$change_vst_median = change_group$vst %>% median()
  entry$neutral_vst_median = neutral_group$vst %>% median()
  
  if("cr_l2fc_50" %in% names(entry)) {
    entry$change_cr_l2fc_50_mean = change_group$cr_l2fc_50 %>% mean()
    entry$neutral_cr_l2fc_50_mean = neutral_group$cr_l2fc_50 %>% mean()
    entry$change_cr_l2fc_50_median = change_group$cr_l2fc_50 %>% median()
    entry$neutral_cr_l2fc_50_median = neutral_group$cr_l2fc_50 %>% median()
  }
  if(!is.null(prefix)){
    entry_names=names(entry)
    names(entry) = str_replace(entry_names,"change_",prefix)
  }
  return(entry)
}

get_wilcox_expression_groups_vst = function(cohort_gene_expression_long,target_gene_lst,gene_id_col="gene_id",
                                            attr_col="call",attr_value="gain",attr_value_contrast="0",prefix=NULL) {
  df = data.frame()
  target_gene_lst = target_gene_lst[target_gene_lst %in% cohort_gene_expression_long[,gene_id_col]]
  
  for(target_gene in target_gene_lst) {
    entry=c()
    
    target_expression = cohort_gene_expression_long[cohort_gene_expression_long[,gene_id_col]==target_gene,]
    if(nrow(target_expression)==0) { next() }
    test_group = target_expression[target_expression[,attr_col]==attr_value,]
    
    if(!is.null(attr_value_contrast)) {
      contrast_group = target_expression[target_expression[,attr_col]==attr_value_contrast,]
    } else {
      contrast_group = target_expression[target_expression[,attr_col]!=attr_value,]
    }
    
    
    if(nrow(test_group)==0|nrow(contrast_group)==0) { next() } 
    
    entry = wilcox_test_groups_vst_cna(test_group,contrast_group,prefix=prefix)
    entry[gene_id_col] = target_gene
    
    df = rbind(df,entry)
    
  } 
  
  if(nrow(df)==0) {
    print("ERROR: no df, check target gene list ")
    return()
  }
  
  if(is.null(prefix)){
    prefix="change_"
  }  
  
  df$change_vst_mean = df[,paste0(prefix,"vst_mean")]
  df$change_vst_median = df[,paste0(prefix,"vst_median")]
  
  df = df  %>% 
    mutate( l2fc_mean=log2(change_vst_mean/neutral_vst_mean),
            l2fc_median=log2(change_vst_median/neutral_vst_median)) 
  
  if(prefix!="change_") {
    #remove tmp
    df = df %>% dplyr::select(-change_vst_mean, -change_vst_median)
  }
  
  df = df %>% arrange(-abs(l2fc_mean))
  return(df)
}


prepare_expression_df_deseq2 = function(expression_overview_wide,patient_lst,filter_read_counts = 10) {
  
  ## subset to patients
  expression_overview_wide = expression_overview_wide[names(expression_overview_wide) %in% c("ensembl_id",patient_lst)]
  
  ## some genes appear to be duplicated but these are the PAR_Y's => remove
  expression_overview_wide = expression_overview_wide[!grepl("PAR_Y",expression_overview_wide$ensembl_id),]
  expression_overview_wide$ensembl_id_noversion =remove_version_from_id(expression_overview_wide$ensembl_id)
  duplicated_genes = expression_overview_wide[duplicated(expression_overview_wide$ensembl_id_noversion),c("ensembl_id_noversion")]
  if(length(duplicated_genes)>1) {
    expression_overview_wide[expression_overview_wide$ensembl_id_noversion %in% duplicated_genes,] 
    print("WARNING duplicated genes")
  } else {
    expression_overview_wide$ensembl_id = expression_overview_wide$ensembl_id_noversion
    expression_overview_wide=expression_overview_wide %>% dplyr::select(-ensembl_id_noversion)
  }
  
  rownames(expression_overview_wide) = expression_overview_wide$ensembl_id
  #renove NAs if exist
  expression_overview_wide[complete.cases(expression_overview_wide),]
  #remove ensembl id 
  expression_overview_wide=expression_overview_wide[-1]
  
  #remove lowly expressed
  expression_overview_wide = expression_overview_wide[rowMeans(expression_overview_wide)>=filter_read_counts,]
  
  return(expression_overview_wide)
  
}

subset_vst = function(genes_vst, sign_genes, gene_display_cols=c("ensembl_id","gene_name")) {
  
  #prevent duplicates
  sign_genes = unique(sign_genes[,gene_display_cols])
  
  sign_genes_vsd = genes_vst[ rownames(genes_vst) %in% sign_genes$ensembl_id, ] %>% as.data.frame()
  sign_genes_vsd$ensembl_id = rownames(sign_genes_vsd)
  sign_genes_vsd = sign_genes_vsd %>% left_join(sign_genes[,gene_display_cols])
  rownames(sign_genes_vsd) = sign_genes_vsd$gene_name
  #sign_genes_vsd %>% dplyr::filter(is.na(gene_name))
  #sign_genes %>% filter(ensembl_id %in% filter(sign_genes_vsd,is.na(gene_name))$ensembl_id)
  sign_genes_vsd = sign_genes_vsd %>% dplyr::select(-all_of(gene_display_cols))
  return(sign_genes_vsd)
}


plot_distances = function(sign_genes_vsd,heatmap_anno,title="",custom_colors=NULL) {
  sampleDists <- dist(t(sign_genes_vsd))
  sampleDistMatrix <- as.matrix( sampleDists )
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  p=pheatmap(sampleDistMatrix,
             annotation_col= heatmap_anno,
             annotation_colors = custom_colors,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors,
             main = title)
  return(p)
}

plot_heatmap = function(input_vsd,heatmap_anno,title="",cluster_rows_on=T,cluster_cols_on=T,custom_colors=NULL,annotation_rows=NULL,custom_color_ramp_palette="RdYlBu",custom_breaks=NA) {
  
  p = pheatmap(mat=input_vsd,
               annotation_col = heatmap_anno,
               annotation_row = annotation_rows,
               annotation_colors = custom_colors,
               cluster_rows = cluster_rows_on,cluster_cols = cluster_cols_on,
               main=title,
               fontsize_row = 8,
               na_col = "darkgrey",
              breaks = custom_breaks,
                  color = colorRampPalette(rev(brewer.pal(n = 7, name =custom_color_ramp_palette)))(100))
  return(p)
}


get_ora_df= function(sign_genes_direction,universe_gene_lst,orgDb,ontology="BP",skip_directional=F,pval=0.05) {
  
  ora_go_df=data.frame()
  if(skip_directional==F){
    for(direction in c("up","down")){
      sign_genes_direction_subset = sign_genes_direction[sign_genes_direction$direction==direction,]
      if(nrow(sign_genes_direction_subset)==0) { next() }
      ora_go = enrichGO( sign_genes_direction_subset$ensembl_id, 
                         universe=universe_gene_lst,
                         OrgDb=orgDb,ont=ontology,
                         keyType = "ENSEMBL",pvalueCutoff = pval, pAdjustMethod = "BH", qvalueCutoff = 0.2)
      
      if(is_null(ora_go)) { next() }
      ora_go_symbols = setReadable(ora_go, OrgDb = orgDb)
      ora_go_subset=ora_go_symbols %>% clusterProfiler::simplify() %>%  as.data.frame() 
      if(nrow(ora_go_subset)==0) { next() }
      ora_go_subset$direction=direction
      ora_go_df=rbind(ora_go_df,ora_go_subset)
    }
  }
  #all
  ora_go_all = enrichGO( sign_genes_direction$ensembl_id, 
                         universe=universe_gene_lst,
                         OrgDb=orgDb,ont=ontology,
                         keyType = "ENSEMBL",pvalueCutoff = pval, pAdjustMethod = "BH", qvalueCutoff = 0.2)
  
  if(!is.null(ora_go_all)) {
    ora_go_symbols = setReadable(ora_go_all, OrgDb = orgDb)
    ora_go_all_df=ora_go_symbols %>% clusterProfiler::simplify() %>%  as.data.frame() 
    if(nrow(ora_go_all_df)>0) {
      ora_go_all_df$direction="both"
      ora_go_df = rbind(ora_go_df,ora_go_all_df)
    }
  }
  
  return(ora_go_df)
}


get_gsea_df = function(sign_genes,orgDb,ontology="BP",value_col="log2FoldChange") {
  #make sorted gene list 
  
  gene_list = sign_genes[,value_col]
  names(gene_list) = sign_genes[,c("ensembl_id")]
  gene_list = sort(gene_list,decreasing = T)
  
  gsea=gseGO(gene_list,
             OrgDb=orgDb,ont=ontology,
             keyType = "ENSEMBL",pvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  gsea_symbols = setReadable(gsea, OrgDb = orgDb)
  gsea_df = gsea_symbols %>% simplify() %>% as.data.frame()
  return(gsea_df)
  
  
  #NB:     #up only and down only does not work wants full genes lists
  
}


get_gse_reactome_df = function(sign_genes,orgDb,value_col="log2FoldChange") {
  ensembl_to_entrez = bitr(sign_genes[,c("ensembl_id")], fromType="ENSEMBL", toType="ENTREZID",OrgDb = orgDb) %>% dplyr::rename(ensembl_id=ENSEMBL)
  sign_genes = sign_genes %>% left_join(ensembl_to_entrez)
  
  sign_genes = sign_genes %>% filter(!is.na(ENTREZID))
  gene_list = sign_genes[,value_col]
  names(gene_list) = sign_genes[,c("ENTREZID")]
  gene_list = sort(gene_list,decreasing = T)
  
  
  reactome=gsePathway(gene_list,organism = "human",
                      #OrgDb=orgDb,
                      pvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  if(is.null(reactome)) {return()}
  reactome_symbols = setReadable(reactome, OrgDb = orgDb)
  reactome_df = reactome_symbols  %>% as.data.frame()
  
  return(reactome_df)
  
}

get_enrich_reactome_df = function(sign_genes,orgDb,universe_gene_lst) {
  ensembl_to_entrez = bitr(sign_genes[,c("ensembl_id")], fromType="ENSEMBL", toType="ENTREZID",OrgDb = orgDb) %>% dplyr::rename(ensembl_id=ENSEMBL)
  sign_genes = sign_genes %>% left_join(ensembl_to_entrez)
  
  sign_genes = sign_genes %>% filter(!is.na(ENTREZID))
  gene_list = sign_genes$ENTREZID
  
  reactome=enrichPathway(gene_list,organism = "human",
                         universe = universe_gene_lst,
                         pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 5)
  
  if(is.null(reactome)) {return()}
  reactome_symbols = setReadable(reactome, OrgDb = orgDb)
  reactome_df = reactome_symbols  %>% as.data.frame()
  
  return(reactome_df)
  
}


get_enrich_wp_df = function(sign_genes,orgDb,universe_gene_lst) {
  ensembl_to_entrez = bitr(sign_genes[,c("ensembl_id")], fromType="ENSEMBL", toType="ENTREZID",OrgDb = orgDb) %>% dplyr::rename(ensembl_id=ENSEMBL)
  sign_genes = sign_genes %>% left_join(ensembl_to_entrez)
  
  sign_genes = sign_genes %>% filter(!is.na(ENTREZID))
  gene_list = sign_genes$ENTREZID
  
  wp=enrichWP(gene_list,organism = "Homo sapiens" ,
              universe = universe_gene_lst,
              pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.2,minGSSize = 5)
  
  if(is.null(wp)) {return()}
  wp_symbols = setReadable(wp, OrgDb = orgDb)
  wp_df = wp_symbols  %>% as.data.frame()
  
  return(wp_df)
  
}
