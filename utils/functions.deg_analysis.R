prepare_expression_df_deseq2 = function(expression_overview_wide,patient_lst,filter_read_counts=10) {
  
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

plot_heatmap = function(input_vsd,heatmap_anno,title="",cluster_rows_on=T,cluster_cols_on=T,
                        custom_colors=NULL,annotation_rows=NULL,fontsize=8,annotation_legend_bool=T,
                        custom_color_ramp_palette="RdYlBu",custom_breaks=NA,gaps_col=NULL,gaps_row=NULL,
                        border_color="grey60",na_color="darkgrey",silent=FALSE) {
  
  p = pheatmap(mat=input_vsd,
               annotation_col = heatmap_anno,
               annotation_row = annotation_rows,
               annotation_colors = custom_colors,
               cluster_rows = cluster_rows_on,cluster_cols = cluster_cols_on,
               main=title,
               fontsize_row = fontsize,
               na_col = na_color,
               drop_levels = TRUE,
               annotation_legend=annotation_legend_bool,
               color=colorRampPalette(rev(brewer.pal(n = 7, name =custom_color_ramp_palette)))(100),
               breaks=custom_breaks,
              gaps_col = gaps_col,gaps_row = gaps_row,
              border_color = border_color, silent=silent,)
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
