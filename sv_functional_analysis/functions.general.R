uq_patients = function(df) {
  return( df %>% dplyr::select(patient_id) %>% unique() %>% nrow())
}


uq_gene_names = function(df) {
  return( df %>% dplyr::select(gene_name) %>% unique() %>% nrow())
}

uq_gene_ids = function(df) {
  return( df %>% dplyr::select(gene_id) %>% unique() %>% nrow())
}

remove_version_from_id = function(identifier) {
  return(sub("\\..*","",identifier))
}

rbind_no_colmatch = function (df1,df2)  {
  if(nrow(df1)==0) {
    return(df2)
  }
  cols_1 = names(df1)
  cols_2 = names(df2)
  col_diff_1 = cols_2[!cols_2 %in% cols_1]
  col_diff_2 = cols_1[!cols_1 %in% cols_2]
  df1[,col_diff_1]=NA
  df2[,col_diff_2]=NA
  
  join_df=rbind(df1,df2)
  return(join_df)
}


cnt_rows = function(df, attr_name=NULL, group_cols=c("patient_id")) {
  cnt_table = df %>% group_by(across(all_of(group_cols))) %>% summarise(svs_cnt = n(),.groups="keep") %>% ungroup()
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c(group_cols,attr_name)
  }
  return(cnt_table)
}


cnt_patients = function(df, group_cols=c("gene_id","gene_name")) {
  cnt_table = df %>% 
    group_by(across(all_of(group_cols))) %>%
    summarize(
      patient_cnt = length(unique(patient_id)),
      patient_lst =toString(sort(unique(patient_id))),
    ) %>%
    arrange(-patient_cnt)
  
  return(cnt_table)
}


print_fraction = function(df, filter_criteria,cnt_subject=c("gene_id")) {
  subset = df %>% dplyr::filter_(filter_criteria)
  
  
  df_cnt = df[,cnt_subject] %>% unique() %>% as.data.frame() %>% nrow()
  subset_cnt = subset[,cnt_subject] %>% unique()  %>% as.data.frame() %>% nrow()
  
  return(paste0(filter_criteria,": ", subset_cnt, " (",round(subset_cnt/df_cnt,2),")"))
  
}
add_missing_patients_to_overview = function(overview_df,cohort){
  missing_patients = cohort %>% dplyr::filter(!patient_id %in% overview_df$patient_id)
  if(nrow(missing_patients)==0) {return(overview_df)}
  missing_cols = names(overview_df)[!names(overview_df) %in% names(missing_patients)]
  missing_patients[,missing_cols]=NA
  overview_df = rbind(overview_df,missing_patients[,names(overview_df) ])
  return(overview_df)
}

## External resources
get_cosmic_genes = function(cosmic_path_input=NULL) {
  if(!is.null(cosmic_path_input)){
    cosmic_path=cosmic_path_input
  }
  cosmic = read.table(cosmic_path,header=T,sep="\t")
  
  cosmic_genes = strsplit(paste(cosmic$Gene.Symbol,cosmic$Synonyms,sep=",",collapse = ","),",")[[1]]
  cosmic_genes = remove_version_from_id(cosmic_genes)
  
  oncogenes = strsplit(paste( dplyr::filter(cosmic,grepl("oncogene",Role.in.Cancer))$Gene.Symbol,
                              dplyr::filter(cosmic,grepl("oncogene",Role.in.Cancer))$Synonyms,sep=",",collapse = ","),",")[[1]]
  oncogenes = remove_version_from_id(oncogenes)
  
  
  tsg = strsplit(paste( dplyr::filter(cosmic,grepl("TSG",Role.in.Cancer))$Gene.Symbol,
                        dplyr::filter(cosmic,grepl("TSG",Role.in.Cancer))$Synonyms,sep=",",collapse = ","),",")[[1]]
  tsg = remove_version_from_id(tsg)
  
  cosmic_genes = as.data.frame(x = cosmic_genes)
  colnames(cosmic_genes)=c("gene_id")
  cosmic_genes = cosmic_genes %>%  dplyr::mutate(oncogene = gene_id %in% oncogenes, tsg = gene_id %in% tsg)
  return(cosmic_genes)
}                          
get_cancer_genes = function(resources_dir=NULL) {
  if(!is.null(resources_dir)){
  library(stringi)
  map_template_vars=c('${resources_dir}'=resources_dir)
  cosmic_path = stri_replace_all_fixed(cosmic_path_template,names(map_template_vars), map_template_vars,vectorize=F)
  oncokb_path = stri_replace_all_fixed(oncokb_path_template,names(map_template_vars), map_template_vars,vectorize=F)
  grobner_recurrent_path = stri_replace_all_fixed(grobner_recurrent_path_template,names(map_template_vars), map_template_vars,vectorize=F)
  }
  
  cosmic_genes=get_cosmic_genes(cosmic_path)
  cosmic_genes$source="cosmic"
  
  cancer_genes = cosmic_genes #gene id, oncogene, tsg
  
  oncokb = read.table(oncokb_path,header=T,sep="\t")
  oncokb = oncokb[,c("gene_name","oncogene","tsg")] %>% 
    mutate(oncogene=ifelse(oncogene=="Yes",T,F), tsg=ifelse(tsg=="Yes",T,F),source="oncokb") %>% dplyr::rename(gene_id=gene_name)
  ## TODO: ensembl tx to gene name, annotate databases that it occurs in?
  # names(oncokb)
  
  cancer_genes=rbind(cancer_genes,oncokb)
  
  grobner_recurrent = read.table(grobner_recurrent_path,header=T,sep="\t") 
  grobner_recurrent=grobner_recurrent %>% dplyr::mutate(gene_id = gene, source="grobner")
  
  grobner_onco = grobner_recurrent[grobner_recurrent$alteration_type=="Amplification",]
  grobner_onco$oncogene=T
  grobner_onco$tsg=NA 
  
  
  grobner_tsg = grobner_recurrent[grobner_recurrent$alteration_type=="Deletion" | 
                                    grobner_recurrent$alteration_type=="Gene-disrupting structural variant" ,]
  grobner_tsg$oncogene=NA
  grobner_tsg$tsg=T
  
  cancer_genes=rbind(cancer_genes, grobner_onco[,names(cancer_genes)], grobner_tsg[,names(cancer_genes)])
  
  cancer_genes = cancer_genes %>% filter(gene_id!="")
  return(cancer_genes)
}

get_gene_onco_tsg_consensus = function(cancer_genes) {
  gene_onco_tsg_consensus = cancer_genes %>% 
    group_by(gene_id) %>% 
    summarize(
      db_lst=toString(sort(unique(source))),
      db_cnt=length(unique(source)),
      onco=sum(oncogene==T),
      tsg=sum(tsg==T)
    )
  
  gene_onco_tsg_consensus[is.na(gene_onco_tsg_consensus)]=0
  gene_onco_tsg_consensus$onco_or_tsg = c("onco","tsg")[max.col(gene_onco_tsg_consensus[,c("onco","tsg")], ties.method = "first")]
  
  gene_onco_tsg_consensus = gene_onco_tsg_consensus %>% mutate(onco_or_tsg = ifelse(onco==0&tsg==0,"unknown",ifelse(onco==tsg,"onco/tsg",onco_or_tsg)))
  
  gene_onco_tsg_consensus$gene_name=gene_onco_tsg_consensus$gene_id
  
  return(gene_onco_tsg_consensus)
}

get_df_setdiff = function(df1,df2,columns_to_compare,key_id) {
  compare_df1 = df1[,columns_to_compare]
  compare_df1=data.frame(lapply(compare_df1, as.character), stringsAsFactors=FALSE)
  
  compare_df2_csr=df2[,columns_to_compare]
  compare_df2_csr=data.frame(lapply(compare_df2_csr, as.character), stringsAsFactors=FALSE)
  
  diff_current_df = dplyr::setdiff(compare_df1,compare_df2_csr) %>% dplyr::rename_with(function(x){paste0("cur_",x)})
  diff_newest_df = dplyr::setdiff(compare_df2_csr,compare_df1) %>% dplyr::rename_with(function(x){paste0("new_",x)})
  
  merged_diff = diff_current_df %>%  merge(diff_newest_df,by.x=paste0("cur_",key_id), by.y=paste0("new_",key_id))
  
  
  #"cur_patient_id"="new_patient_id"))
  return(merged_diff)
}


get_genome_windows = function(target_width=1e6,fixed_size=F) {
  contig_lengths_path = paste0(resources_dir,"contig_lengths.tsv")
  contig_lengths = read.table(contig_lengths_path,header=T,sep="\t")
  contig_lengths$end=contig_lengths$length
  contig_lengths$start=0
  contigs = GRanges(contig_lengths)
#  contig_tiles = tile(contigs,width=target_width) %>% unlist() 
#default 

  if(fixed_size) {
    bins = slidingWindows(contigs,width=target_width,step=target_width) %>% unlist()
  } else {
    bins = tile(contigs,width=target_width) %>% unlist()
  }
  
  return(bins)
}


get_most_variable_rows = function(df,ntop=500) {
  # select the ntop rows by variance
  if (is.data.frame(df)) {
    df = df %>% keep(is.numeric)  
    df = as.matrix(df)
  }
  
  rv <- rowVars(df)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(df)))]
  #amonst the top are XIST and TSIXT => subset to autosomes 
  
  top_rows = df[select,]  
  top_rows = as.data.frame(top_rows)
  return(top_rows)
}

get_gene_properties_df = function(gtf_path) {
  gtf <- rtracklayer::import(gtf_path)
  gene_properties = gtf[gtf$type=="gene"]
  names(gene_properties) = gene_properties$gene_id
  gene_properties_df = as.data.frame(gene_properties)
  gene_properties_df$to_coordinate = paste0(gene_properties_df$seqnames,":",gene_properties_df$start,"-",gene_properties_df$end)
  gene_properties_df=gene_properties_df %>% dplyr::rename(gene_width=width)
  gene_properties_df$ensembl_id=remove_version_from_id(gene_properties_df$gene_id)
  return(gene_properties_df)
}

get_exon_properties_df = function(gtf_path){
  gtf <- rtracklayer::import(gtf_path)
  exon_properties = gtf[gtf$type=="exon"&!is.na(gtf$exon_id) & (gtf$gene_type=="protein_coding" | grepl("RNA",gtf$gene_type))]
  exon_properties$exon_row = paste0("exon_",1:length(exon_properties))
  
  exon_properties_df = as.data.frame(exon_properties)
  exon_properties_df$to_coordinate = paste0(exon_properties_df$seqnames,":",exon_properties_df$start,"-",exon_properties_df$end)
  exon_properties_df$ensembl_id=remove_version_from_id(exon_properties_df$gene_id)
  return(exon_properties_df)
}


get_genes_to_cytobands = function(gene_properties,chromosome_bands) {
  genes_to_cytobands = get_reciprocal_overlap_pairs(gene_properties,chromosome_bands,svtype_matching = F,reciprocal_overlap = 0)
  genes_to_cytobands$gene_id=genes_to_cytobands$set1
  genes_to_cytobands$cytoband=genes_to_cytobands$set2
  #be inclusive to boundary genes
  #genes_to_cytobands %>% filter(overlap_set1_set2<1)
  
  gene_to_cytoband_mapping = genes_to_cytobands %>% group_by(gene_id) %>% 
    summarize(cytoband=toString(unique(sort(cytoband)))) %>% ungroup() %>% as.data.frame()
  return(gene_to_cytoband_mapping)
}


## also in functions.svs.df => can I remove here?
cnt_svs = function(svs_df, attr_name=NULL, group_cols=c("patient_id")) {
  cnt_table = svs_df %>% group_by(across(all_of(group_cols))) %>% summarise(svs_cnt = n(),.groups="keep") %>% ungroup()
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c(group_cols,attr_name)
  }
  return(cnt_table)
}


get_replication_timing_df  = function(replication_timing_path) {
  replication_timing_df = read.table(replication_timing_path,header=T,sep="\t",comment.char = "")
  
  replication_timing_df = replication_timing_df %>% dplyr::rename(seqnames=CHR,start=POS,rt_anno=switch) %>% dplyr::mutate(end=start+49999) 
  replication_timing_df$rt_id= paste0("rt_",1:nrow(replication_timing_df))
  replication_timing_df$to_coordinate = paste0(replication_timing_df$seqnames,":",replication_timing_df$start,"-",replication_timing_df$end)
  
  return(replication_timing_df)
}



add_missing_patients_to_overview = function(overview_df,cohort){
  missing_patients = cohort %>% dplyr::filter(!patient_id %in% overview_df$patient_id)
  if(nrow(missing_patients)==0) {return(overview_df)}
  missing_cols = names(overview_df)[!names(overview_df) %in% names(missing_patients)]
  missing_patients[,missing_cols]=NA
  missing_patients=unique(missing_patients[,names(overview_df) ])
  overview_df = rbind(overview_df,missing_patients)
  return(overview_df)
}


get_recurrent_gene_centric = function(gene_centric_df,grouping_cols=c("gene_name")) {
  recurrent_df = gene_centric_df %>% group_by(across(all_of(grouping_cols))) %>% 
    summarise(patients= toString(unique(patient_id)), patient_cnt=length(unique(patient_id))) %>%
    arrange(-patient_cnt)
  
  recurrent_df = recurrent_df %>% as.data.frame()
  return(recurrent_df)
}



annotate_gene_group = function(df,cancer_genes,wilms_genes_of_interest=NULL) {
  if(!is.null(wilms_genes_of_interest)) {
    df = df %>% mutate(
      gene_group = ifelse(gene_name %in% filter(cancer_genes,oncogene==T)$gene_id & gene_name %in% filter(cancer_genes,tsg==T)$gene_id,"onco/tsg",
                          ifelse(gene_name %in% filter(cancer_genes,oncogene==T)$gene_id,"oncogene",
                                 ifelse(gene_name %in% filter(cancer_genes,tsg==T)$gene_id,"tsg",       
                                        ifelse(gene_name %in% wilms_genes_of_interest$gene_name,"wilms_other",
                                               ifelse(gene_name %in% cancer_genes$gene_id,"cancer_other",NA))))))
  } else {
    df = df %>% mutate(
      gene_group = ifelse(gene_name %in% filter(cancer_genes,oncogene==T)$gene_id & gene_name %in% filter(cancer_genes,tsg==T)$gene_id,"onco/tsg",
                          ifelse(gene_name %in% filter(cancer_genes,oncogene==T)$gene_id,"oncogene",
                                 ifelse(gene_name %in% filter(cancer_genes,tsg==T)$gene_id,"tsg",       
                                        ifelse(gene_name %in% cancer_genes$gene_id,"cancer_other",NA)))))
    
  }
  
  return(df)
}


#normalize for plotting 

df_normalize_for_plotting = function(df) {
  rownames(df) = df$patient_id
  df = df %>% dplyr::select(-patient_id) %>% as.data.frame()
  df = df/rowSums(df)
  df$patient_id = rownames(df)
  return(df)
}

# Total (within) sum of squares for a cluster
cluster_within_sum_square = function(matrix, group) {
  sum(aggregate(matrix, by=list(group), function(x) sum(scale(x,scale=FALSE)^2))[, -1])
}


matrix_edit_cell = function(mat,filter_df,rowname_attr="gene_name_display",col_attr="patient_id",value=NA) {
  target_lst = filter_df[,rowname_attr] %>% unique()
  target_lst = target_lst[target_lst %in% rownames(mat)]
  
  for(row in target_lst) {
    target_col_lst = filter_df[filter_df[,rowname_attr]==row,col_attr]
    target_col_lst = target_col_lst[target_col_lst %in% colnames(mat)]
    
    mat[row,target_col_lst]=value
  }
  return(mat)
}

get_patient_gene_summary = function(df,attr_name=NULL) {
  summary = df %>% group_by(patient_id) %>%
    summarize(genes = toString(unique(sort(paste0(gene_name," (",alteration," ",ifelse(!is.na(call),call,""),")")))))
  
  if(!is.null(attr_name)) {
    summary[,attr_name]=summary[,c("genes")]
    summary = summary %>% dplyr::select(-genes)
  }
  return(summary)
}


