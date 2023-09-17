uq_patients = function(df) {
  if("patient_label" %in% names(df)) {
    return( df %>% dplyr::select(patient_label) %>% unique() %>% nrow())
  }
  if("patient_id" %in% names(df)) {
    return( df %>% dplyr::select(patient_id) %>% unique() %>% nrow())
  }
}


uq_gene_names = function(df) {
  return( df %>% dplyr::select(gene_name) %>% unique() %>% nrow())
}

uq_gene_ids = function(df) {
  return( df %>% dplyr::select(gene_id) %>% unique() %>% nrow())
}

uq_svs = function(merged_svs) {
  return(length(unique(merged_svs$patient_sv_merged)))
}


cnt_uq = function(df, uq_id="patient_sv_merged", attr_name=NULL, group_cols=c("patient_label")) {
  cnt_table = df %>% group_by(across(all_of(group_cols))) %>% summarise(svs_cnt = length(unique(!! sym(uq_id))),.groups="keep") %>% ungroup()
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c(group_cols,attr_name)
  }
  return(cnt_table)
}


lst_str = function(charstring){
  return(toString(sort(unique(charstring))))
}


cnt_str = function(charstring){
  return(length(unique(charstring)))
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


subset_maintain_label_order = function(bins,patient_label_lst) {
  subset_bins = bins %>% filter(patient_label %in% patient_label_lst)
  subset_bins$patient_label = factor(subset_bins$patient_label,level=unique(bins$patient_label))
  return(subset_bins)
}

get_csr_subject_info = function(csr_export_path) {
  csr_export=read.table(csr_export_path,sep="\t",header=T)
  
  csr_subject_info = csr_export %>% dplyr::mutate(
    patient_id=Subject.Id,
    diagnosis_id=Diagnosis.Id,
    patient_label=str_replace(diagnosis_id,"PMCDN","M"),
    diagnosis_date=X04..Date.of.diagnosis,
    sex=X03..Sex,
    age_at_first_diagnosis = X01..Age.at.first.diagnosis,
    birth_date=X02..Date.of.birth,
    death_date=X06..Date.of.death.)
  
  csr_subject_info = csr_subject_info %>% dplyr::select(patient_id,patient_label,diagnosis_id,sex,age_at_first_diagnosis,death_date,birth_date,diagnosis_date) %>% unique()
  return(csr_subject_info)
}


resize_gr_distance = function(gr,distance,flag_as_distance=TRUE) {
  if(flag_as_distance==TRUE){
    distance=distance*2
  }
  return(GenomicRanges::resize(gr,width=width(gr)+distance,fix = "center"))
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
add_missing_patients_to_overview = function(overview_df,cohort,flag_patient_label = F){
  if(flag_patient_label==T){
    missing_patients = cohort %>% dplyr::filter(!patient_label %in% overview_df$patient_label)
  } else {
    missing_patients = cohort %>% dplyr::filter(!patient_id %in% overview_df$patient_id)
  }
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



get_genes_of_interest_collection = function() {
  ## just for quick loading in betwen documents without creating an export
  #lot of work to do to clean this up!
  # eg special document with genes for cancer types?
  #         flag_gene_of_interest = gene_name %in% genes_of_interest_collection$gene_name,
  #         flag_pmc_panel = gene_name %in% filter(genes_of_interest_collection,grepl("pmc_",db_lst))$gene_name) 

  cancer_genes = get_cancer_genes(resources_dir)
  gene_onco_tsg_consensus = get_gene_onco_tsg_consensus(cancer_genes)
  
  
  fusion_list_diagnostics = read.table(paste0(resources_dir,"InclusionList_fusion_genes_Homo_sapiens_GRCh38_BioMart_v1.8.txt"),sep="\t",header=T)
  
  gene_panel_diagnostics = read.table(paste0(resources_dir,"diagnostic_somatic_2.5.lst"),sep="\t",header=F)
  gene_list_diagnostics = c(fusion_list_diagnostics$Gene.name,gene_panel_diagnostics$V1) %>% unique()
  #gene_list_diagnostics[!gene_list_diagnostics %in% cancer_related$gene_name]
  
  #Wnt pathway
  library(msigdbr)
  msigdb_c2  = msigdbr::msigdbr(species="Homo sapiens",category="C2")
  wnt_pathway = msigdb_c2 %>% filter(gs_id=="M39669")
  wnt_pathway_genes = wnt_pathway %>% dplyr::rename(gene_name=gene_symbol, ensembl_id = ensembl_gene)
  #wnt_pathway_genes = wnt_pathway_genes %>% left_join(gene_properties_df[,gene_cols])
  
  neuroblastoma_genes_of_interest = c("FOXR1","PTPRD", "CSMD1", 
                                      "ODZ2","ODZ3", "ODZ4","TENM2","TENM3","TENM4",
                                      "FBXO8","CEP44",
                                      "MYCN","ALK","ATRX","TERT",
                                      "TFAP2B", "MAP7", "PTPRH","SLC18A1",
                                      "KIF1B","PLEKHG5","UBE4B", "CHD5", "CADM1", "ATM",
                                      "RASGRP3", "SMARCE1", "EBF1", "SDHA", "HMGA2",
                                      "NRAS","HRAS","KRAS","BRAF","PTPN11", "NF1") #molenaar elevelt ras genes
  
  
  wilms_genes_of_interest = c("LIN28A","LIN28B","MIRLET7A")
  
  actionable_mutations_villani = read.table(paste0(resources_dir,"actionable_mutations.Villani_2022.txt"),sep = "\t",header=T)
  
  #construct
  genes_of_interest_collection = gene_onco_tsg_consensus %>% dplyr::select(gene_name,db_lst,onco_or_tsg)
  
  genes_of_interest_collection = rbind_no_colmatch(genes_of_interest_collection,
                                                   neuroblastoma_genes_of_interest %>% as.data.frame() %>% dplyr::rename(gene_name=".") %>% dplyr::mutate(db_lst="neuroblastoma_manual"))
  genes_of_interest_collection = rbind_no_colmatch(genes_of_interest_collection,
                                                   wilms_genes_of_interest %>% as.data.frame() %>% dplyr::rename(gene_name=".") %>% dplyr::mutate(db_lst="wilms_manual"))
  
  genes_of_interest_collection = rbind_no_colmatch(genes_of_interest_collection,
                                                   wnt_pathway_genes %>% dplyr::rename(db_lst=gs_name) %>% dplyr::select(gene_name,db_lst) )
  
  genes_of_interest_collection = rbind_no_colmatch(genes_of_interest_collection,
                                                   fusion_list_diagnostics %>% dplyr::mutate(gene_name=Gene.name,db_lst="pmc_fusions") %>% dplyr::select(gene_name,db_lst) )
  
  
  genes_of_interest_collection = rbind_no_colmatch(genes_of_interest_collection,
                                                   gene_panel_diagnostics %>% dplyr::mutate(gene_name=V1,db_lst="pmc_somatic_panel") %>% dplyr::select(gene_name,db_lst) )
  
  genes_of_interest_collection = rbind_no_colmatch(genes_of_interest_collection,
                                                   actionable_mutations_villani %>% dplyr::mutate(db_lst="actionable_mutations_villani") %>% dplyr::select(gene_name,db_lst) )
  
  
  #use both:
  #fusion_list_diagnostics %>% filter(!Gene.name %in% gene_panel_diagnostics$V1)
  return(genes_of_interest_collection)
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


get_chr_arms = function(chromosome_bands_df,split_giestain=T,autosomes= c(paste("chr",1:22,sep=""),"chrX","chrY")) {
  ## note; 1 based coordinates??? 
  
  chromosome_bands_df = chromosome_bands_df %>% 
    dplyr::filter(seqnames %in% autosomes) %>% 
    dplyr::mutate(chr_arm = ifelse(grepl("p",cytoband),paste0(seqnames,"p"),paste0(seqnames,"q")))
  
  if(split_giestain==T){
    chromosome_bands_df = chromosome_bands_df %>% 
      dplyr::mutate(chr_arm = ifelse(grepl("gneg",gieStain)|grepl("gpos",gieStain),chr_arm,paste0(chr_arm,"_",gieStain))) 
  }
  chr_arms = chromosome_bands_df %>% dplyr::group_by(seqnames,chr_arm) %>% 
    dplyr::summarize(start = min(start)+1,
                     end = max(end)) %>% as.data.frame() %>% ungroup()
  
  return(chr_arms %>% as.data.frame())
}

get_chromosomes = function(df) {
  ## NB: cannot remove the acen/gvar regions from chromosomes because it spans the centromere
  #start is always 1 but depends on giestain used or not so this is more flexible
  
  if("chr_arm" %in% names(df)) {
    chr_arms_df = df
  } else if ("cytoband" %in% names(df)) {
    chr_arms_df = get_chr_arms(df)
  } else {
    return("Please input chromosome bands or chr arms dataframe")
  }
  
  chromosomes_df = chr_arms_df %>%  group_by(seqnames) %>% summarize(start = min(start),end=max(end),.groups="drop")
  chromosomes = GRanges(chromosomes_df)
  names(chromosomes)=chromosomes_df$seqnames
  chromosomes_df$chrom_width=chromosomes_df$end-chromosomes_df$start+1 #to make equal with chr arm width
  return(chromosomes_df %>% as.data.frame())
}

get_genome_bins = function(target_width=1e6,fixed_size=F,contig_lengths_path=NULL) {
  if(is.null(contig_lengths_path)) {
    contig_lengths_path = paste0(resources_dir,"contig_lengths.tsv")
  }
  contig_lengths = read.table(contig_lengths_path,header=T,sep="\t")
  contig_lengths$end=contig_lengths$length
  contig_lengths$start=0
  contigs = GRanges(contig_lengths)

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

get_recurrent_gene_centric = function(gene_centric_df,grouping_cols=c("gene_name"),flag_patient_label = F) {
  
  
  if(flag_patient_label) {
    recurrent_df = gene_centric_df %>% group_by(across(all_of(grouping_cols))) %>% 
      summarise(patients= toString(unique(patient_label)), patient_cnt=length(unique(patient_label))) %>%
      arrange(-patient_cnt)
  } else {
    recurrent_df = gene_centric_df %>% group_by(across(all_of(grouping_cols))) %>% 
      summarise(patients= toString(unique(patient_id)), patient_cnt=length(unique(patient_id))) %>%
      arrange(-patient_cnt)
  }
  
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

add_missing_bins = function(cnt_by_bins,cohort,genome_bins_df) {
  
  all_bins_patients = tidyr::crossing(patient_label=cohort$patient_label, bin_id=genome_bins_df$bin_id) 
  cnt_by_bins_incl_missing = cnt_by_bins
  missing_bins=anti_join(all_bins_patients,cnt_by_bins)
  missing_bins = missing_bins %>% 
    left_join(genome_bins_df[,c("bin_id","seqnames")]) %>%
    mutate(cnt=0)
  
  cnt_by_bins_incl_missing = rbind_no_colmatch(missing_bins,cnt_by_bins_incl_missing)
  cnt_by_bins_incl_missing$bin_id = factor(cnt_by_bins_incl_missing$bin_id,
                                           levels=unique(genome_bins_df[gtools::mixedorder(genome_bins_df$bin_id),]$bin_id))
  cnt_by_bins_incl_missing$seqnames = factor(cnt_by_bins_incl_missing$seqnames,levels=unique(genome_bins_df[ gtools::mixedorder(genome_bins_df$seqnames),c("seqnames")]))
  
  
  return(cnt_by_bins_incl_missing)
}


BinomPvalsNrAlts <- function(cnts, nalts, ngenes){
  return(sapply(cnts, FUN = function(el){return(binom.test(x = el, n = nalts, p = 1/ngenes, alternative = "greater")$p.value)}))
}

BinomPvalsNrAltsCorLength <- function(df, cntscol, genelengthscol, n, totallength){
  res <- apply(df[,c(cntscol, genelengthscol)], MARGIN = 1, FUN = function(row){
    res2 <- binom.test(x = row[1], n = n, p = row[2]/totallength, alternative = "greater")$p.value
  })
  return(res)
}



get_baseline_correction = function(cohort,automated_baseline_shifts,manual_baseline_shifts){
  names(automated_baseline_shifts)
  baseline_shifts = cohort[,c("patient_id","tumor_id")] %>% left_join(automated_baseline_shifts,by=c("tumor_id"="biomaterial_id","patient_id")) %>% left_join(manual_baseline_shifts) 
  
  baseline_shifts$abs_baseline_correction = as.numeric(baseline_shifts$abs_baseline_correction) #cooerce to NA if other sign
  baseline_shifts = baseline_shifts %>% mutate(baseline_correction = ifelse(!is.na(abs_baseline_correction),abs_baseline_correction, #use manual
                                                                            ifelse(abs(median_neutral_cr)<0.05 & abs(baseline_shift_observed)<0.05,0, #0 if smaller than 0.05
                                                                                   ifelse(is.na(median_neutral_cr),-1*baseline_shift_observed, #if no median neutral then take observed
                                                                                          -1*median_neutral_cr)))) ## TODO this goes wrong often as well. need to check
  baseline_shifts = baseline_shifts %>% mutate(baseline_correction = ifelse(abs(baseline_correction)<0.05,0,baseline_correction))
  
  baseline_shifts$baseline_correction = as.numeric(baseline_shifts$baseline_correction)
  
  return(baseline_shifts)
}


get_gr_coordinate=function(df,attr_name="coordinate") {
  df[,attr_name] = paste0(df$seqnames,":",df$start,"-",df$end)
  return(df)
}
  
get_merged_svs = function(svs_df,sv_merge_cols=c("sv_merged","sv_merged_coordinate"),return_partnered=F) {
  #note, only safe per patient!
  sv_merged_df = svs_df %>% dplyr::rename(tumor_af_bp = tumor_af, normal_af_bp = normal_af) %>%
    group_by(across(all_of(sv_merge_cols))) %>% summarize(
      pass_all = all(FILTER=="PASS"),
      FILTER = toString(unique(sort(FILTER))),
      sv_cnt = length(unique(as.character(sv_name))),
      sv_names = toString(unique(sort(as.character(sv_name)))),
      partners = toString(unique(sort(as.character(partner)))),
      coordinates = toString(unique(sort(coordinate))),
      partner_coordinates = toString(unique(sort(partner_coordinate))),
      
      tools=toString(unique(sort(tool))),
      tools_cnt = length(unique(tool)),
      tumor_af=mean(tumor_af_bp,na.rm=T), normal_af=mean(normal_af_bp,na.rm=T),
      tumor_af_spread=(max(tumor_af_bp,na.rm=T)-min(tumor_af,na.rm=T)), normal_af_spread=(max(normal_af_bp,na.rm=T)-min(normal_af_bp,na.rm=T)), 
      
      svtype=toString(unique(sort(svtype))), 
      svlen=ifelse(!grepl("CTX",svtype),mean(svLen,na.rm=T),NA),
      svlen_spread=ifelse(!is.na(svlen),(max(svLen,na.rm=T)-min(svLen,na.rm = T)),NA),
      .groups="keep") %>% ungroup() 
  
  
  ## show partnered bps
  if(return_partnered==T) {
    print("WARNING only safe per patient!")
    
    partnered = sv_merged_df[,c("sv_merged","partners","sv_names")] %>% 
      left_join(sv_merged_df[,c("sv_merged","partners","sv_names")],
                by=c("partners"="sv_names","sv_names"="partners")) %>% 
      dplyr::rename(partner_sv_merged= sv_merged.y, sv_merged=sv_merged.x)
    
    
    sv_merged_df = sv_merged_df %>% left_join(partnered[,c("sv_merged","partner_sv_merged")],by = "sv_merged")
    
    partner_sv_merged_coord = sv_merged_df[,c("sv_merged","partner_sv_merged","sv_merged_coordinate")] %>%
      left_join(sv_merged_df[,c("sv_merged","partner_sv_merged","sv_merged_coordinate")],
                by=c("partner_sv_merged"="sv_merged","sv_merged"="partner_sv_merged")) %>% 
      dplyr::rename(partner_sv_merged_coordinate= sv_merged_coordinate.y, sv_merged_coordinate=sv_merged_coordinate.x)
    
    sv_merged_df=sv_merged_df %>% left_join(partner_sv_merged_coord[,c("sv_merged","partner_sv_merged","partner_sv_merged_coordinate")],by=c("sv_merged","partner_sv_merged"))
    
  }
  sv_merged_df = sv_merged_df %>% as.data.frame() %>% unique()
  return(sv_merged_df)
}

get_ctx_partner_merged_coordinates = function(svs_df) {
  
  merged_incl_partner_coord_df = data.frame()
  for(pid in unique(svs_df$patient_label)) {
    merged_incl_partner_coord = get_merged_svs(svs_df, return_partnered = T)
    
    merged_incl_partner_coord$patient_label=pid
    merged_incl_partner_coord_df = rbind(merged_incl_partner_coord_df,merged_incl_partner_coord)
  }
  
  merged_incl_partner_coord_df$patient_sv_merged = paste0(merged_incl_partner_coord_df$patient_label,"_",merged_incl_partner_coord_df$sv_merged)
  merged_incl_partner_coord_df$partner_sv_merged = paste0(merged_incl_partner_coord_df$patient_label,"_",merged_incl_partner_coord_df$partner_sv_merged)
  
  merged_incl_partner_coord_df = merged_incl_partner_coord_df[,c("patient_label","patient_sv_merged","partner_sv_merged","partner_sv_merged_coordinate")]
  
  merged_incl_partner_coord_df = unique(merged_incl_partner_coord_df)
  return(merged_incl_partner_coord_df)
  
  
}


make_merged_svs = function(unfiltered_svs_df) {
  merged_svs = unfiltered_svs_df %>%
    #   head() %>% 
    group_by(patient_label,svtype,patient_sv_merged,sv_merged_coordinate,flag_in_filtered_svs) %>% summarize(
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
  
  #we want partner coordinates > do merging per patient 
  merged_incl_partner_coord_df = get_ctx_partner_merged_coordinates(filter(unfiltered_svs_df,svtype=="CTX"))
  
  merged_svs = merged_svs %>% left_join(merged_incl_partner_coord_df,by=c("patient_label","patient_sv_merged"))
  return(merged_svs)
}

load_cohort_svs = function(cohort,map_template_vars,svs_path_template=NULL,flag_add_patient_label=T) {
  if(is.null(svs_path_template)) {
    svs_path_template = svs_union_path_template
  }
  if(is.null(map_template_vars[['${merged_svs_dir}']])) {
    print("merged svs dir undefined, exiting")
    return(NULL)
  }  
  cohort_sv = data.frame()
  for(pid in cohort$patient_id) {
    patient = cohort %>% filter(patient_id==pid)
    map_template_vars_patient = c(map_template_vars,'${patient_basename}'=patient$basename)  
    svs_union_path = stri_replace_all_fixed(svs_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    
    if(length(Sys.glob(svs_union_path))==1) { 
      #Load sv anno
      sv_anno = read.table(svs_union_path,header=T,sep="\t",stringsAsFactors = F)
      if(nrow(sv_anno)==0) { next() }
      cohort_sv = rbind_no_colmatch(cohort_sv,sv_anno)
    } else {
      print(paste0("Missing for patient ",patient$patient_label,patient$basename,": ",svs_union_path))
      next()
    }  
    
  }                      
  
  cohort_sv = cohort_sv %>% left_join(cohort[,c("patient_id","patient_label")])  
  if(flag_add_patient_label) {
    cohort_sv$patient_sv_merged = paste0(cohort_sv$patient_label,"_",cohort_sv$sv_merged)
    cohort_sv$patient_sv_name = paste0(cohort_sv$patient_label,"_",cohort_sv$sv_name)
    cohort_sv$partner_sv_name = paste0(cohort_sv$patient_label,"_",cohort_sv$partner) #for multitool ctx
  }
  return(cohort_sv)
}


annotate_sv_multitool_support = function(svs_df,sv_tumor_af_threshold=0.1) {
  unfiltered_svs_df=svs_df
  ## for all coding analyses filter by tumor AF >0.1 and then reassess multi tool support 
  svs_df = unfiltered_svs_df %>% filter(tumor_af>sv_tumor_af_threshold)
  filtered_svs_multi_tool_support = get_multi_tool_support(svs_df)
  svs_df = svs_df %>% filter(patient_sv_merged %in% filtered_svs_multi_tool_support$patient_sv_merged) 
  unfiltered_svs_df = unfiltered_svs_df %>% mutate(flag_in_filtered_svs=patient_sv_name %in% svs_df$patient_sv_name)
  
  #svs_df %>% nrow() == unfiltered_svs_df %>% filter(flag_in_filtered_svs) %>% nrow()
  return(unfiltered_svs_df)
}


load_cohort_snvs = function(cohort,map_template_vars,analysis_type="somatic",snv_path_template=NULL,flag_add_patient_label=T) {
  if(is.null(snv_path_template)) {
    snv_path_template = snv_output_path_template
  }
  
  if(!"${utils_output_dir}" %in% names(map_template_vars)) {
    print("utils output dir undefined, exiting")
    return(NULL)
  }
  if(!"${analysis_type}" %in% names(map_template_vars) & !is.null(analysis_type)) {
    print(paste0("Setting analysis type ",analysis_type))
    map_template_vars[["${analysis_type}"]]=analysis_type
  }
  
  if(map_template_vars[["${analysis_type}"]]!=analysis_type) {
    print("Template vars and analysis type mismatch")
    print(paste0("map_template_vars: ",map_template_vars[["${analysis_type}"]],"; analysis type:",analysis_type))
  }
  
  
  cohort_snv_df = data.frame()
  for(pid in cohort$patient_id) {
    patient = filter(cohort,patient_id==pid)
    
    if(analysis_type=="somatic") { 
      map_template_vars_patient=c(map_template_vars, '${patient_basename}'=patient$basename)
    } else {
      map_template_vars_patient=c(map_template_vars, '${patient_basename}'=patient$normal_id)
    }
    snv_path = stri_replace_all_fixed(snv_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    
    if(length(Sys.glob(snv_path))==0) { 
      print(paste0("WARNING: missing patient ",patient$patient_id," : ",snv_path))
      next() 
    }
    
    snv_patient = read.table(snv_path,header=T,sep="\t")
    
    version_info=readLines(snv_path,n=1)
    version_info = version_info %>% strsplit(split=" ",fixed=T) %>% unlist()
    vep_version=gsub("[^[:alnum:] ]", "", version_info[1])
    snv_patient$vep_version=vep_version
    
    cohort_snv_df=rbind_no_colmatch(cohort_snv_df,snv_patient) 
    #note we can lose some cols this way, but necessary because of different VEP versions
  }
  
  if(nrow(cohort_snv_df)==0) {
    print("WARNING: empty dataframe!")
    return(cohort_snv_df)
  }
  if(cohort %>% uq_patients() != cohort_snv_df %>% uq_patients()) {
    print("WARNING: patients missing?")
    print(paste0(cohort %>% uq_patients() ," in df: ", cohort_snv_df %>% uq_patients()))
  }
  
  # TODO future if CDS mapping is back then 
  #mutate(display = ifelse(!is.na(cds_mapping),cds_mapping,snv_id))
  cohort_snv_df = cohort_snv_df %>% dplyr::rename(gene_name=SYMBOL,ensembl_id=Gene) %>% mutate(display=snv_id)
  
  if(flag_add_patient_label) {
    cohort_snv_df = cohort_snv_df %>% left_join(cohort[,c("patient_id","patient_label")])
    cohort_snv_df$patient_snv_id=paste0(cohort_snv_df$patient_label,"_",cohort_snv_df$snv_id)
  }
  return(cohort_snv_df)
}


load_cohort_cn_segments = function(cohort,map_template_vars,segments_path_template=NULL,baseline_shifts=NULL,flag_add_patient_label=T,flag_silence_messages=F) {
#  cn_segments_path_template= paste0("${cna_data_dir}/${biomaterial_id}*_${sequencing_strategy}*${cna_seg_file_ext}")
  
  if(is.null(segments_path_template)) {
    segments_path_template=cn_segments_path_template
  }
  
  if(!"${cna_data_dir}" %in% names(map_template_vars)) {
    map_template_vars[["${cna_data_dir}"]]=cna_data_dir
  }
  
  if(!"${sequencing_strategy}" %in% names(map_template_vars)) {
    map_template_vars[["${sequencing_strategy}"]]="WGS"
  }
  if(!"${cna_seg_file_ext}" %in% names(map_template_vars)) {
    map_template_vars[["${cna_seg_file_ext}"]]=cna_seg_file_ext
  }
  
  segments_df = data.frame()
  
  for(pid in cohort$patient_label) {
    patient = cohort %>% filter(patient_label == pid)
    
    map_template_vars_patient=c(map_template_vars, '${biomaterial_id}'=patient$tumor_id)
    segments_df_path = stri_replace_all_fixed(segments_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
    
    segments_df_file=Sys.glob(segments_df_path)
    if(length(segments_df_file)!=1){ 
      if(flag_silence_messages==F) {
        print("WARNING: input unclear")
        print(segments_df_path)
      }
      next()
    }
    
    contig_lengths = get_contig_lengths(segments_df_file)
    
    if(sum(contig_lengths$length) != expected_autosomal_length) {
      if(flag_silence_messages==F) {
        print("WARNING: autosomal length different than expected, is it really hg38?")
      }
    }
    
    modeled_seg= read_modeled_seg(segments_df_file)
    #todo refactor get_gr_coordinate()
    modeled_seg$from_coordinate = paste0(modeled_seg$seqnames,":",modeled_seg$start,"-",modeled_seg$end)
    
    
    #default
    #modeled_seg_cols=c("cna_id","width","cr_l2fc_50","maf_50","cr_l2fc_10","maf_10","cr_l2fc_90","maf_90")  
    modeled_seg_cols=c(modeled_seg_cols,"from_coordinate")
    modeled_seg_cols = modeled_seg_cols[modeled_seg_cols %in% names(modeled_seg)]
    
    modeled_seg$patient_label=patient$patient_label 
    modeled_seg$tumor_id=patient$tumor_id
    
    segments_df = rbind(segments_df,modeled_seg)
  }
  if(nrow(segments_df)==0) {
    return(data.frame())
  }
  if(!is.null(baseline_shifts)) {
    segments_df = segments_df %>% left_join(baseline_shifts)
    segments_df[is.na(segments_df$baseline_correction),c("baseline_correction")]=0
  } else {
    segments_df = segments_df %>% left_join(cohort[,c("patient_id","tumor_id")])
    segments_df$baseline_correction=0  
  }
  
  segments_df$cr_l2fc_50 = segments_df$cr_l2fc_50+segments_df$baseline_correction
  
  
  segments_df = segments_df %>% call_cna()
  segments_df = segments_df %>% get_gr_coordinate("cna_coordinate")
  
  ## TODO refactor with label
  segments_df$patient_cna_id=paste0(segments_df$patient_id,"_",segments_df$cna_id)
  
  if(flag_add_patient_label) {
    segments_df = segments_df %>% left_join(cohort[,c("patient_id","patient_label")])
  }
  return(segments_df)
}



plot_pca = function(matrix,annotation,  identifier_col = "patient_id",color_by=NULL,shape_by=NULL,title="",pc_x=1,pc_y=2,flag_return_empty=F) {
  pca = prcomp(matrix)
  percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
  
  if(FALSE) { 
    print(summary(pca))
    
    rotations = as.data.frame(pca$rotation)
    rotations$ensembl_id = rownames(rotations)
    rotations = rotations %>% left_join(gene_properties_df,by="ensembl_id") %>% arrange(-abs(PC1)) 
    #write.table(rotations,pca_rotations_path,sep="\t",quote=F,row.names = F,col.names = T)
    
    top_genes_pc1= rotations %>% arrange(-abs(PC1)) 
    top_genes_pc2= rotations %>% arrange(-abs(PC2)) 
  }
  
  pcaData = as.data.frame(pca$x)
  pcaData[,identifier_col] =rownames(pcaData)
  pcaData=merge(pcaData, annotation,by=identifier_col)

    p=ggplot(data=pcaData, aes_string(x = paste0("PC",pc_x), y = paste0("PC",pc_y)))   +
    xlab(paste0("PC",pc_x,": ", percentVar[pc_x], "% variance"))+ylab(paste0("PC",pc_y,": ", percentVar[pc_y], "% variance")) +
    ggtitle(title) +
    theme_bw() 
    
    if(flag_return_empty==F) {
      p=p+geom_point(aes_string(color=color_by,shape=shape_by),size=4)  
    }
  return(p)
}


plot_umap = function(matrix,umap_max=5,annotation,identifier_col="patient_id",color_by=NULL,shape_by=NULL,title="",umap_min=1,flag_return_empty=F) {
  pca = prcomp(matrix)
  umap_max = min(umap_max,length(as.data.frame(pca$x)))
  
  umat = umap(as.data.frame(pca$x)[,umap_min:umap_max]) 
  umat = as.data.frame(umat$layout[,1:2])
  colnames(umat) = c("dim1", "dim2")
  umat[,identifier_col] = rownames(umat)
  umat = merge(umat, annotation, by=identifier_col)
  
  p=ggplot(data=umat, aes_string(x = "dim1", y = "dim2")) +
    theme_bw() +
    ggtitle(paste0("UMAP of PCA ",umap_min,"-",umap_max," of ",title))
  
  if(flag_return_empty==F) {
    p=p+geom_point(aes_string(color=color_by,shape=shape_by),size=4)  
  }
  
  return(p)
}

plot_binned_cn = function(cn_by_bins,color_cancer_type=T) {
  p = ggplot(
    cn_by_bins %>% filter(!grepl("acen|stalk|gvar",bin_id)),aes(x=bin_id,y=patient_label)) +
    geom_tile(
      aes(x=bin_id,fill=call,y=patient_label),size=1,width=1, height=0.9) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=annotation_colors_call,na.value = "white") +
    theme(axis.text.x=element_blank())+
    scale_x_discrete(breaks = NULL) +  
    facet_grid(cols=vars(seqnames),space = "free",scales = "free")
  
  if(color_cancer_type) {
    p=p+scale_color_manual(values= c(annotation_colors$cancer_type)) 
  }
  return(p)
}

plot_cnpare_grouping = function(p,type="EUC") {
  
  #for(k in c(3,6,8)) {}  vars(eval(!!paste0("EUC_k",k))) doesnt work yet
  if(type=="EUC") {
    print(p +  facet_grid(rows=vars(EUC_k3),cols=vars(seqnames),space = "free",scales = "free"))
    print(p +  facet_grid(rows=vars(EUC_k6),cols=vars(seqnames),space = "free",scales = "free"))
    print(p +  facet_grid(rows=vars(EUC_k8),cols=vars(seqnames),space = "free",scales = "free"))
  }
  if(type=="MAN") {
    print(p +  facet_grid(rows=vars(MAN_k3),cols=vars(seqnames),space = "free",scales = "free") )
    print(p +  facet_grid(rows=vars(MAN_k6),cols=vars(seqnames),space = "free",scales = "free") )
    print(p +  facet_grid(rows=vars(MAN_k8),cols=vars(seqnames),space = "free",scales = "free") )
  }
  
  if(type=="COR") {
    print(p +  facet_grid(rows=vars(COR_k4),cols=vars(seqnames),space = "free",scales = "free") )
  }
  if(type=="COS") {
    print(p +  facet_grid(rows=vars(COS_k5),cols=vars(seqnames),space = "free",scales = "free") )
  }
  
}

plot_binned_cn_sv_ctx_taf = function(cn_sv_bp_taf_by_bins,taf_sv_by_bins_long,cn_ctx_by_bins,color_cancer_type=T,plot_patient_label_tumor=F) {
  
  if(plot_patient_label_tumor) {
    cn_sv_bp_taf_by_bins = cn_sv_bp_taf_by_bins %>% select(-patient_label) %>% dplyr::rename(patient_label=patient_label_tumor)
    taf_sv_by_bins_long = taf_sv_by_bins_long %>% select(-patient_label) %>% dplyr::rename(patient_label=patient_label_tumor)
    cn_ctx_by_bins = cn_ctx_by_bins %>% select(-patient_label) %>% dplyr::rename(patient_label=patient_label_tumor)
  }
  
  p = ggplot(
    cn_sv_bp_taf_by_bins %>% filter(!grepl("acen|stalk|gvar",bin_id)),aes(x=bin_id,y=patient_label)) +
    geom_tile(
      aes(x=bin_id,fill=call,y=patient_label),size=1,width=1, height=0.9) + 
    geom_tile(
      data=filter(taf_sv_by_bins_long,cnt>0 & !grepl("acen|stalk|gvar",bin_id))  , 
      aes(x=bin_id,y=patient_label), size=1,width=2,height=0.3,fill="darkgrey",alpha=0.5) +
    geom_tile(
      data=filter(cn_ctx_by_bins,sv_overlap & !grepl("acen|stalk|gvar",bin_id)), 
      aes(x=bin_id,y=patient_label),size=1,width=2, height=0.5,fill="black") + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=annotation_colors_call,na.value = "white") +
    theme(axis.text.x=element_blank())+
    facet_grid(cols=vars(seqnames),space = "free",scales = "free")+
    scale_x_discrete(breaks = NULL) 
  
  if(color_cancer_type) {
    p=p+scale_color_manual(values= c(annotation_colors$cancer_type)) 
  }
  return(p)
}

count_gr_start_end_by_bins = function(gr,genome_bins,gr_df=data.frame(),gr_id_col="patient_sv_merged",group_cols="patient_label",flag_apply_ordering=T) {
  overlaps = get_reciprocal_overlap_pairs_start_end(gr,genome_bins,reciprocal_overlap = 0,svtype_matching = F)
  overlaps = overlaps %>% dplyr::rename(!!gr_id_col := set1,bin_id=set2)
  
  overlaps = overlaps %>% left_join(gr_df[,c(group_cols,gr_id_col,"seqnames")],by=gr_id_col)
  
  by_bins = overlaps %>% 
    group_by(across(c(group_cols,"seqnames","bin_id"))) %>% 
    dplyr::summarize(cnt=length(unique(!! sym(gr_id_col)))) %>% as.data.frame()
  
  
  if(flag_apply_ordering){
    by_bins$bin_id = factor(by_bins$bin_id,levels=unique(genome_bins[gtools::mixedorder(genome_bins$bin_id)]$bin_id))
    by_bins$seqnames = factor(by_bins$seqnames,levels=unique(chr_arms_df[ gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")]))
  }
  return(by_bins)
}

## CN/SV integration scripts ----


get_split_range_by_breakpoint = function(svs,ranges,ranges_df,patients_svs_df,svs_id_col="patient_sv_merged",ranges_id_col="chr_arm"){
  
  ## Make window start -- SV -- end 
  ## can generalize to other ranges like chr arms if you switch out 'seqnames' and 'chr arms'
  
  #assumptions  
  #sv start < sv end
  #range start < chr arm end 
  
  # before =  range start to SV start
  # after = SV end to range end
  
  svs_coordinates_df = patients_svs_df %>% dplyr::select(c(all_of(svs_id_col),"seqnames","start","end")) %>% 
    pivot_longer(cols=c("start","end"),names_to="sv_breakpoint_orientation",values_to="sv_coordinate")
  
  ranges_coordinates_df = ranges_df %>% dplyr::select(all_of(unique(c(ranges_id_col,"seqnames","start","end")))) %>% 
    pivot_longer(cols=c("start","end"),names_to="orientation",values_to="chr_arm_coordinate")
  
  svs_range_overlaps = get_reciprocal_overlap_pairs_start_end(svs,ranges,reciprocal_overlap = 0,svtype_matching = F)
  
  svs_range_overlaps = svs_range_overlaps %>% dplyr::rename(!!svs_id_col:=set1,!!ranges_id_col:=set2) %>% 
    left_join(svs_coordinates_df[,c(svs_id_col,"sv_coordinate","sv_breakpoint_orientation")],by=c(svs_id_col,"sv_breakpoint_orientation")) %>%
    left_join(ranges_coordinates_df,by=c(ranges_id_col,"sv_breakpoint_orientation"="orientation"))
  
  #to before/after ranges
  windows_df = svs_range_overlaps %>% 
    mutate(
      start = ifelse(sv_breakpoint_orientation=="start",chr_arm_coordinate,sv_coordinate),
      end = ifelse(sv_breakpoint_orientation=="start",sv_coordinate,chr_arm_coordinate),
      window_orientation = ifelse(sv_breakpoint_orientation=="start","before","after"))
  
  windows_df$window_id=paste0(windows_df[,svs_id_col],"_",windows_df$window_orientation)
  windows_df = windows_df %>% get_gr_coordinate("window_coordinate")
  windows_df$window_width = windows_df$end-windows_df$start
  return(windows_df)
}

annotate_stable_size_pass_window = function(df, unbalanced_ctx_window_width_min=10e6,min_window_width_both=T,stable_windows_both=T){
  #require size minimum for cn change or not, set to 0 then always true
  #smallest chrom is 46 mb seems reasonable to require both to be 10 mb minimum
  if(is.null(unbalanced_ctx_window_width_min)) {
    unbalanced_ctx_window_width_min=0
  }
  if(min_window_width_both) {
    df = df %>% mutate(flag_stable_size_pass = 
                         frac_cr_stable_before & frac_cr_stable_after & 
                         window_width_before > unbalanced_ctx_window_width_min & window_width_after > unbalanced_ctx_window_width_min)
    
  } else if(stable_windows_both & min_window_width_both==F & unbalanced_ctx_window_width_min>0) { 
    #one or the other
    df = df %>% mutate(flag_stable_size_pass = 
                         frac_cr_stable_before & frac_cr_stable_after &  (
                           window_width_before > unbalanced_ctx_window_width_min | window_width_after > unbalanced_ctx_window_width_min))
  } else if(stable_windows_both==F & min_window_width_both==F) {
    df = df %>% mutate(flag_stable_size_pass = 
                         (frac_cr_stable_before & window_width_before > unbalanced_ctx_window_width_min) | 
                         (frac_cr_stable_after & window_width_after > unbalanced_ctx_window_width_min))
    
  } else {
    print("this shouldnt happen")
  }
  
  
  return(df)
}
annotate_cn_change_window = function(df, cr_l2fc_threshold_sv_change=0.2,apply_cr_value_minimum=F){
  #2023-05-05 update: flag_cn_change_window can be T/F or NA depending on flag_stable_size_pass, if not pass then always NA
  #feature_value_cn_change= balanced if no change and stable size pass
  df = df %>% 
    dplyr::mutate(feature_value_before = paste0(frac_call_before, " (cr:",round(frac_cr_l2fc_50_before,2),", ",round(window_width_before/1e6,2)," Mbp, stable:", round(frac_seg_covered_cr_stable_before,2),")"),
                  feature_value_after = paste0( frac_call_after, " (cr:",round(frac_cr_l2fc_50_after,2),", ",round(window_width_after/1e6,2)," Mbp, stable:", round(frac_seg_covered_cr_stable_after,2),")"))
  
  if(!"flag_stable_size_pass" %in% names(df)) {
    print("need flag_stable_size_pass")
    return(df)
  }
  
  #call change sufficient?
  if(apply_cr_value_minimum==F) {
    df = df %>% 
      dplyr::mutate(flag_cn_change_window = ifelse(flag_stable_size_pass, 
                                                   ( frac_call_before != frac_call_after | abs(frac_cr_l2fc_50_before-frac_cr_l2fc_50_after) > cr_l2fc_threshold_sv_change), NA))
    
  } else if(apply_cr_value_minimum==T) {
    df = df %>% 
      dplyr::mutate(flag_cn_change_window = ifelse(flag_stable_size_pass, 
                                                   (abs(frac_cr_l2fc_50_before-frac_cr_l2fc_50_after) > cr_l2fc_threshold_sv_change),NA))
  }
  ## not sure if feature_value_cn_change makes sense like this?
  ## not always the right approach see example of gain-0 with bp in q arm  M158AAB_merged_81a3d6f825f4aa0187fd053ab43efc78
  if(FALSE) {
    df = df %>% 
      mutate(feature_value_cn_change = ifelse(flag_cn_change_window, 
                                              ifelse(grepl("p",chr_arm), 
                                                     paste0(feature_value_before," -- ",frac_call_after), 
                                                     paste0(frac_call_before," -- ",feature_value_after)),NA))
  }
  df = df %>% mutate(feature_value_cn_change = ifelse(!is.na(flag_cn_change_window),
                                                      ifelse(flag_cn_change_window,paste0(frac_call_before,"-",frac_call_after), "balanced"),NA))
  return(df)
}


make_window_around_sv_bp = function(svs_gr,sv_bp_flank_window_size=1e6,svs_id_col="patient_sv_merged") {
  #flank before start breakpoint
  flank_svs_1mb_start=svs_gr
  strand(flank_svs_1mb_start)="*" #remove strand for uniform flanks
  flank_svs_1mb_start = flank(flank_svs_1mb_start,sv_bp_flank_window_size,start = TRUE)
  all( svs_gr %>% start() -  flank_svs_1mb_start %>% start() == 
         sv_bp_flank_window_size) & all( svs_gr  %>% start() -  flank_svs_1mb_start %>% end() == 1 )
  
  svs_1mb_start_coord = flank_svs_1mb_start %>% as.data.frame() %>% get_gr_coordinate("window_coordinate")
  svs_1mb_start_coord$window_orientation="before"
  svs_1mb_start_coord$window_id = paste0(svs_1mb_start_coord[,svs_id_col],"_before")
  
  #flank after end breakpoint
  flank_svs_1mb_end=svs_gr
  strand(flank_svs_1mb_end)="*" #remove strand for uniform flanks
  flank_svs_1mb_end = flank(flank_svs_1mb_end,sv_bp_flank_window_size,start = FALSE)
  (all(svs_gr  %>% end() -  flank_svs_1mb_end %>% start() == -1) & 
      all(svs_gr %>% end() -  flank_svs_1mb_end %>% end() == -sv_bp_flank_window_size))
  
  svs_1mb_end_coord = flank_svs_1mb_end %>% as.data.frame() %>% get_gr_coordinate("window_coordinate")
  svs_1mb_end_coord$window_orientation="after"
  svs_1mb_end_coord$window_id = paste0(svs_1mb_end_coord[,svs_id_col],"_after")
  
  svs_1mb_windows = rbind(svs_1mb_start_coord,svs_1mb_end_coord)
  svs_1mb_windows$window_width = svs_1mb_windows$width
  return(svs_1mb_windows)
}


assess_cn_state_window_before_after = function(segments,segments_df,ctx_windows_chrom_df,windows_df_cols=c("window_id","window_coordinate","window_width"),svs_id_col="patient_sv_merged",patient_id_col="patient_label",threshold_window_frac_cn_state=0.7){
  
  ctx_windows_chrom_gr = GRanges(ctx_windows_chrom_df$window_coordinate)
  names(ctx_windows_chrom_gr) = ctx_windows_chrom_df$window_id
  mcols(ctx_windows_chrom_gr) = ctx_windows_chrom_df[,c(patient_id_col,windows_df_cols)]
  
  window_chrom_cn_fractions = get_chrom_cn_fractions(segments,segments_df,
                                                     ranges=ctx_windows_chrom_gr,
                                                     ranges_df=ctx_windows_chrom_df,
                                                     cohort,return_wide = F,
                                                     ranges_id_col="window_id",ranges_width_col="window_width",relative_cr = F)
  
  window_chrom_cn_fractions_highest = window_chrom_cn_fractions[
    window_chrom_cn_fractions$frac_covered == ave(window_chrom_cn_fractions$frac_covered, window_chrom_cn_fractions$window_id,FUN=max),]
  #assign window_frac_stable based on threshold stable state
  window_chrom_cn_fractions_highest = window_chrom_cn_fractions_highest %>% 
    dplyr::mutate(cr_stable = seg_covered_cr_stable > threshold_window_frac_cn_state) 
  
  window_chrom_cn_fractions_highest = window_chrom_cn_fractions_highest %>% ungroup() %>%  select(all_of(window_cn_fractions_cols))
  #note that this mean cr is for the call state fraction not overall.
  
  #make wide and link up
  ctx_windows_chrom_wide = ctx_windows_chrom_df %>% 
    select(all_of(c(svs_id_col,"window_orientation",windows_df_cols))) %>% 
    pivot_wider(names_from="window_orientation",values_from=all_of(windows_df_cols))
  
  
  
  ctx_chrom_unfiltered = ctx_windows_chrom_wide %>%
    left_join(window_chrom_cn_fractions_highest %>%
                dplyr::rename_with(.cols=all_of(window_cn_fractions_cols_pivot),.fn=function(x){paste0("frac_",x)})  %>%
                dplyr::rename_with(.cols=all_of(c("window_id",paste0("frac_",window_cn_fractions_cols_pivot))),.fn=function(x){paste0(x,"_before")})) %>% 
    left_join(window_chrom_cn_fractions_highest %>%
                dplyr::rename_with(.cols=all_of(window_cn_fractions_cols_pivot),.fn=function(x){paste0("frac_",x)})  %>%
                dplyr::rename_with(.cols=all_of(c("window_id",paste0("frac_",window_cn_fractions_cols_pivot))),.fn=function(x){paste0(x,"_after")})) %>% 
    as.data.frame() 
  
  return(ctx_chrom_unfiltered)
}
if(FALSE) {
get_windows_between_svs = function(svs_df,svs_id_col="patient_sv_merged",chr_id_col="chrom"){
  if(nrow(svs_df)==0) { return(data.frame()) }
  sv_combinations = expand.grid(svs_df[,svs_id_col],svs_df[,svs_id_col])
  
  local_coord_cols = c("start","end")
  sv_combinations = sv_combinations %>% dplyr::rename(self=Var1, other=Var2) %>% 
    left_join(svs_df[,c(svs_id_col,chr_id_col,local_coord_cols)],by=c("self"=svs_id_col)) %>%
    dplyr::rename_with(.cols=local_coord_cols,.fn=function(x){paste0("self_",x)}) %>%
    left_join(svs_df[,c(svs_id_col,chr_id_col,local_coord_cols)],by=c("other"=svs_id_col)) %>% 
    dplyr::rename_with(.cols=local_coord_cols,.fn=function(x){paste0("other_",x)}) 
  
  ## what svs to compare?
  #select svs same chrom and 
  #remove overlapping from the comparison list 
  #because also == the combo with self is removed as well.
  #SELF start > OTHER start & SELF start < OTHER end: 
  #look if SELF start is within OTHER body, remove within and back overlap
  #SELF start < OTHER start & SELF end > OTHER start
  #look if OTHER start is within the self SV body, remove spanning and front overlap. 
  
  sv_combinations = sv_combinations %>%
    filter(chrom.x==chrom.y) %>% dplyr::rename(seqnames=chrom.x) %>% dplyr::select(-chrom.y) %>%
    filter( !(self_start >= other_start & self_start <= other_end) & 
              !(self_start <= other_start & self_end >= other_start) 
    )  
  
  if(nrow(sv_combinations)==0) { return(data.frame()) }
  
  #always have 1 window between SELF and OTHER:
  #  if SELF end < OTHER start #self is before other and window is right/after self 
  #SELF end -- OTHER start 
  #if SELF start > OTHER end #self is after other and window is left/before self 
  #OTHER end --- SELF start 
  #else: situations not thought of?.. 
  
  sv_combinations = sv_combinations %>% rowwise() %>% 
    mutate(
      start = ifelse(self_end<other_start, self_end, ifelse(self_start>other_end,other_end,NA)),
      end = ifelse(self_end<other_start, other_start, ifelse(self_start>other_end,self_start,NA)),
      window_orientation = ifelse(self_end<other_start, "after", ifelse(self_start>other_end,"before",NA)),
      window_id = paste(sort(c(self,other)),collapse = "_")) %>% as.data.frame() 
  
  #window ids show that window from A to B is same as from B to A 
  # sv_combinations %>% select(seqnames,start,end,window_id) %>% unique()
  
  sv_combinations = sv_combinations %>% get_gr_coordinate("window_coordinate")
  sv_combinations$window_width = sv_combinations$end-sv_combinations$start
  
  return(sv_combinations)
}
get_windows_before_after = function(windows_df,window_analysis_cols=c("window_orientation","call","cr_stable"), return_source_df=T,cr_l2fc_threshold_sv_change=0.2) {
  if(return_source_df & "ctx_windows_before_after_change" %in% names(windows_df)) {
    print("Warning, column ctx_windows_before_after_change was already present in windows df and now removed")
    windows_df = windows_df %>% dplyr::select(-ctx_windows_before_after_change)
  }
  window_analysis_cols_pivot=window_analysis_cols[!window_analysis_cols %in% c("window_orientation")]
  
  #before after and see if they differ
  cna_windows_before_after = windows_df %>% filter(window_orientation %in% c("before","after"))
  if(nrow(cna_windows_before_after)>0) {
    cna_windows_before_after = cna_windows_before_after %>% 
      select(all_of(c("patient_sv_merged",window_analysis_cols))) %>% 
      pivot_wider(names_from="window_orientation",values_from=all_of(window_analysis_cols_pivot)) %>%
      mutate(ctx_windows_before_after_change = ( (call_before %in% c("gain","loss") | call_after %in% c("gain","loss")) ) &
             ( (call_before!=call_after) | ( abs(cr_l2fc_50_before-cr_l2fc_50_after) > cr_l2fc_threshold_sv_change) ) )
    
  } else {
    print("Before/After not found in window_orientation")
  }
  
  if(return_source_df) {
    windows_df = windows_df %>% left_join(cna_windows_before_after)
    return(windows_df)
  } else {
    return(cna_windows_before_after)
  }
}
}


get_windows_gr_from_ctx_properties = function(ctx_window_properties, region="chr_arm",property_cols=c("")) {
  # verified this version used for 2023-07 versions
  #before or after windows depending on p/q arm
  ctx_window_properties= ctx_window_properties %>% as.data.frame() %>% ungroup()
  if(!region %in% c("chrom","chr_arm")) {
    print("make sure region is chrom/chr_arm")
    return()
  }
  call_before=paste0("window_",region,"_frac_call_before")
  call_after=paste0("window_",region,"_frac_call_after")
  coordinate_before=paste0("window_",region,"_coordinate_before")
  coordinate_after=paste0("window_",region,"_coordinate_after")
  id_before=paste0("window_",region,"_id_before")
  id_after=paste0("window_",region,"_id_after")
  
  windows_ctx_before = ctx_window_properties %>% filter(chr_arm_orientation=="p")
  windows_ctx_after = ctx_window_properties %>% filter(chr_arm_orientation=="q")
  
  # check if none missing and non both
  missing_ctx = ctx_window_properties %>% filter(
    !patient_sv_merged %in% windows_ctx_before$patient_sv_merged &
      !patient_sv_merged %in% windows_ctx_after$patient_sv_merged)
  
  if(missing_ctx %>% nrow() > 0) {
    print("No windows selected for: ")
    print(toString(missing_ctx$patient_sv_merged))
    return()
  }
  
  windows_ctx_before_gr = GRanges(windows_ctx_before[,coordinate_before])
  names(windows_ctx_before_gr) = windows_ctx_before[,id_before]
  mcols(windows_ctx_before_gr) = windows_ctx_before[,c("chr_arm_orientation","patient_label")]
  
  windows_ctx_after_gr = GRanges(windows_ctx_after[,coordinate_after])
  names(windows_ctx_after_gr) = windows_ctx_after[,id_after]
  mcols(windows_ctx_after_gr) = windows_ctx_after[,c("chr_arm_orientation","patient_label")]
  
  windows_ctx_gr = c(windows_ctx_before_gr,windows_ctx_after_gr)
  windows_ctx_gr$window_id=names(windows_ctx_gr)
  windows_ctx_df = mcols(windows_ctx_gr) %>% as.data.frame()
  windows_ctx_df$region = region
  
  properties_p = ctx_window_properties %>% filter(chr_arm_orientation=="p")
  properties_p = properties_p %>% dplyr::mutate(call := !!sym(call_before), window_id = !!sym(id_before)) %>% select(window_id,call,all_of(property_cols))
  
  properties_q = ctx_window_properties %>% filter(chr_arm_orientation=="q")
  properties_q = properties_q %>% dplyr::mutate(call := !!sym(call_after), window_id = !!sym(id_after)) %>% select(window_id,call,all_of(property_cols))
  
  windows_ctx_df = windows_ctx_df %>% left_join(rbind(properties_p,properties_q))
  mcols(windows_ctx_gr) =windows_ctx_df
  
  return(windows_ctx_gr)
}

get_segments_windows_overlaps = function(segments,segments_df,windows_gr,modeled_seg_cols=c("patient_cna_id","width","cr_l2fc_50","maf_50"),windows_df=NULL,window_cols=NULL) {
  #get the segs overlapping windows needing only window gr made with get_windows_gr_from_ctx_properties()
  segments_windows_overlaps = get_reciprocal_overlap_pairs(segments, windows_gr, svtype_matching = F,reciprocal_overlap = 0,ignore_strand=T)
  segments_windows_overlaps = segments_windows_overlaps %>% dplyr::rename(patient_cna_id=set1,window_id=set2)
  
  segments_windows_overlaps = segments_windows_overlaps %>% 
    left_join(segments_df[,unique(c("patient_label",modeled_seg_cols))],by=c("patient_cna_id"))
  
  if(!is.null(window_cols) & !is.null(windows_df)) {
    segments_windows_overlaps = segments_windows_overlaps %>% 
      left_join(windows_df[,unique(c("patient_label",window_cols))],by="window_id",relationship = "many-to-many") 
    segments_windows_overlaps = segments_windows_overlaps %>% dplyr::rename(patient_label=patient_label.x)
  } else {
    #not needed bc of patient label in window id
    segments_windows_overlaps = segments_windows_overlaps %>% rowwise() %>% dplyr::mutate(patient_label.y=unlist(str_split(window_id,"_"))[1]) 
  }
  
  segments_windows_overlaps = segments_windows_overlaps %>% filter(patient_label==patient_label.y) %>% dplyr::select(-patient_label.y)
  segments_windows_overlaps = segments_windows_overlaps %>% filter(overlap_set1_set2>0.9)
  return(segments_windows_overlaps %>% as.data.frame())
}

# Functions from recurrent svs.Rmd also includes for complex svs ----

plot_single_svs = function(target_svs,plot_title=NULL,start_coord=NULL,end_coord=NULL,display_intrasv_labels=F,display_ctx_labels=T,plot_order=NULL,collapse_svtype=F) {
  
  annotation_colors$svtype = c(`DUP`="red3",`DEL`="steelblue3",`INV`= "mediumpurple2",`CTX`="purple3")
  
  svs_coord = GRanges(target_svs$sv_merged_coordinate)
  mcols(svs_coord)= target_svs
  svs_coord_gr = svs_coord
  svs_coord = svs_coord %>% as.data.frame() #%>% left_join(patient_id_to_labels)
  svs_coord = svs_coord %>% mutate(height=ifelse(svtype=="CTX",0.6,0.5))
  
  #maintain order
  if(is.null(plot_order)){
    svs_coord$patient_sv_merged = factor(svs_coord$patient_sv_merged,levels=target_svs$patient_sv_merged)
  } else {
    svs_coord$patient_sv_merged = factor(svs_coord$patient_sv_merged,levels=plot_order)
    if(collapse_svtype) {
      svs_coord$svtype = factor(svs_coord$svtype,levels=plot_order)
    }
  }
  
  if(is.null(start_coord)) {
    start_coord= min(filter(svs_coord)$start)-100
  }
  if(is.null(end_coord)) {
    end_coord= max(filter(svs_coord)$end)+100
  }
  p = ggplot(svs_coord,
             aes(xmin=(start/1e6),xmax=(end/1e6),ymin=-(height/2),ymax=(height/2),fill=svtype)) +  
    geom_hline(yintercept=0,color="darkgrey") + 
    geom_rect(data=svs_coord %>% filter(svtype!="CTX"), show.legend = T,alpha=0.3,color="black",size=0.1) + 
    geom_rect(data=svs_coord %>% filter(svtype=="CTX"), show.legend = T,alpha=1,color=annotation_colors$svtype["CTX"],size=0.5) + 
    scale_fill_manual(values=c(annotation_colors$svtype),na.value = "white") + 
    theme_classic() + theme(axis.text.y = element_blank()) + xlab("Coordinate in Mb") + ylab("") + 
    ggtitle(plot_title)
  
  if(collapse_svtype) {
    p = p + facet_grid(rows=vars(patient_label,svtype),space = "free_y") 
  } else {
    #each sv its own row
    p = p + facet_grid(rows=vars(patient_label,patient_sv_merged),space = "free_y") 
  }
  
  if(display_ctx_labels) {
    p = p+ geom_label_repel(data=svs_coord %>% filter(svtype=="CTX"), aes(x=start/1e6,y=height,label=paste0(partner_chrom," ",round(tumor_af_mean,2))),size=3, color="white",min.segment.length = 0,nudge_y=0.1,max.overlaps = 30) 
  }
  
  if(display_intrasv_labels) {
    p = p +
      geom_label_repel(data=svs_coord %>% filter(svtype!="CTX"), aes(x=start/1e6,y=(-height/2),label=paste0(svtype," ",svlen_mean_str," ",round(tumor_af_mean,2)),color=svtype),size=3, fill="white",min.segment.length = 0,nudge_y = -0.2,max.overlaps = 30) +
      scale_color_manual(values=c(annotation_colors$svtype),na.value = "white") 
  }
  
  if(!is.null(start_coord)&!is.null(end_coord)&start_coord!="auto"){
    p = p + coord_cartesian(xlim=c(start_coord/1e6,end_coord/1e6))
  }
  
  return(p)  
}
add_start_end_lines = function(p,start_chr_arm,end_chr_arm) {
  p = p +
    geom_vline(xintercept=start_chr_arm/1e6,color="black",size=0.5) +
    geom_vline(xintercept=end_chr_arm/1e6,color="black",size=0.5)
  p=p+theme(panel.spacing = unit(0, "lines"))
  return(p)
}
plot_target_svs_cn_data = function(target_svs,segments_df,cn_around_sv_bp,flag_plot_segval=F,start_coord = NULL,end_coord = NULL,display_intrasv_labels=F,display_ctx_labels=T,collapse_svtype=F,plot_order=NULL,flag_segval_scale_amp=T,display_segval_labels=F,cr_amplification_threshold=3) {
  annotation_colors$svtype = c(`DUP`="red3",`DEL`="steelblue3",`INV`= "mediumpurple2",`CTX`="purple3")
  
  p = plot_single_svs(target_svs,start_coord = start_coord,end_coord = end_coord,display_intrasv_labels = display_intrasv_labels,plot_order = plot_order,display_ctx_labels=display_ctx_labels,collapse_svtype=collapse_svtype)
  p = p + 
    scale_fill_manual(values=c(annotation_colors$svtype,annotation_colors_call),na.value = "white") +
    facet_grid(rows=vars(patient_label,patient_sv_merged,seqnames),space = "free_y")
  
  if(!is.null(cn_around_sv_bp)) {
    #Plot before/after 1mb regions of breakpoints for each SV
    cn_data = cn_around_sv_bp %>% filter(patient_sv_merged %in% target_svs$patient_sv_merged)
    if(cn_data %>% nrow() >0) {
      cn_data$height=0.2
      #cn_data$patient_sv_merged #exists!
      if(!is.null(plot_order)) {
        cn_data$patient_sv_merged = factor(cn_data$patient_sv_merged, levels=plot_order) # as sep. track
      } else {
        cn_data$patient_sv_merged = factor(cn_data$patient_sv_merged,levels = target_svs$patient_sv_merged)
      }
      p = p + geom_rect(data=cn_data, show.legend = T,alpha=1,size=0.1,aes(fill=call),color="black") 
    }
  }
  if(!is.null(segments_df)){
    #plot segments
    cn_data = segments_df[,
                          c("patient_label","tumor_id","patient_cna_id","seqnames","start","end","cr_l2fc_50","call")] %>%
      filter(patient_label  %in% target_svs$patient_label) %>% 
      filter(seqnames %in% target_svs$chrom) %>% #& start >= start_coord & end <= end_coord)  %>%
      unique() 
    
    if(cn_data %>% nrow() >0) {
      cn_data = cn_data %>% call_cna() %>% get_relative_cr_l2fc(cohort) %>% rename_relative_cr_to_normal()
      cn_data[cn_data$cr_l2fc_50>cr_amplification_threshold,"call"]="amp"
      
      cn_data$height=0.2
      cn_data$patient_sv_merged=NULL # overlapping with the SVs
      cn_data$patient_sv_merged = "cn_segments"
      if(!is.null(plot_order) & "cn_segments" %in% plot_order) {
        cn_data$patient_sv_merged = factor(cn_data$patient_sv_merged, levels=plot_order) # as sep. track
      }
      if(collapse_svtype) {
        cn_data$svtype = cn_data$patient_sv_merged
      }
      if(is.null(plot_order) | "cn_segments" %in% plot_order) {
        p = p + geom_rect(data=filter(cn_data,call!="0"), show.legend = T,alpha=1,size=0.1,aes(fill=call)) +
        geom_rect(data=filter(cn_data,call=="0"), show.legend = T,alpha=1,size=0.1,aes(fill=call),color="black")
      }
      
      if(flag_plot_segval) {
        plot_segval = cn_data
        plot_segval[plot_segval$cr_l2fc_50<(-1),"cr_l2fc_50"]=-1
        if(flag_segval_scale_amp) { plot_segval[plot_segval$cr_l2fc_50>cr_amplification_threshold,"cr_l2fc_50"]=cr_amplification_threshold }
        plot_segval[plot_segval$cr_l2fc_50>cr_amplification_threshold,"call"]="amp"
        plot_segval$patient_sv_merged="copy_ratio" # as sep. track
        if(!is.null(plot_order) & "copy_ratio" %in% plot_order) {
          plot_segval$patient_sv_merged = factor(plot_segval$patient_sv_merged, levels=plot_order) # as sep. track
        }
        if(collapse_svtype) {
          plot_segval$svtype = plot_segval$patient_sv_merged
        }
        
        
        p = p + geom_rect(data=plot_segval, aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=cr_l2fc_50,fill=call))
        
        
        if(display_segval_labels) {
          plot_segval$height=plot_segval$cr_l2fc_50
          if(is.null(start_coord) | is.null(end_coord)) {
            target_svs_region = target_svs %>% group_by(seqnames) %>%
              summarize(start_coord = min(start)-100,
                        end_coord = max(end)+100) %>% as.data.frame()
            
            plot_segval = plot_segval %>% left_join(target_svs_region,by="seqnames") %>% filter(start>start_coord,end<end_coord)
          } else {
            plot_segval = filter(plot_segval,start>start_coord,end<end_coord)
          }
          
          p = p +
            geom_label_repel(data=plot_segval, aes(x=start/1e6,y=(cr_l2fc_50),label=paste0("cr:",round(cr_l2fc_50,2))),color="black",size=3, fill="white",min.segment.length = 0,nudge_y = -0.2,max.overlaps = 30) 
        }
        
        
        p = p + facet_grid(rows=vars(patient_label,patient_sv_merged,seqnames),space = "free_y",scale="free_y")
      }
    }
  }
  
  if(collapse_svtype) {
    p = p + facet_grid(rows=vars(patient_label,seqnames,svtype),space = "free_y")
    
  }
  return(p)
}

print_fraction = function(df, filter_criteria,count_rows=T) {
  subset = df %>% dplyr::filter_(filter_criteria)
  if(count_rows) {
    return(paste0(filter_criteria,": ", nrow(subset), " (",round(nrow(subset)/nrow(df),2),")"))
  } 
}


get_min_max_median = function(df,col,attr_prefix=NULL) {
  df = df %>% summarize(min=min(!! sym(col),na.rm = T),
                        median=median(!! sym(col),na.rm = T),
                        mean=mean(!! sym(col),na.rm=T),
                        max=max(!! sym(col),na.rm = T))
  
  if(!is.null(attr_prefix)) {
    df = df %>% dplyr::rename_with(.cols=c("min","max","median","mean"),.fn=function(x){paste0(attr_prefix,x)})
  }
  
  return(df)
}


#todo integrate in other functions
get_range_as_bp_gr = function(range_gr) {
  range_start = range_gr
  end(range_start)=start(range_start)
  range_start$sv_bp_orientation="start"
  
  range_end = range_gr
  start(range_end)=end(range_end)
  range_end$sv_bp_orientation="end"
  
  range_as_bp_gr = c(range_start,range_end)
  return(range_as_bp_gr)
}


gr_midpoint = function(x) {
  start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
  return(x)
}


get_svs_as_bp_gr = function(merged_svs_gr,ctx_midpoint=T) {
  #intra svs split in start and end, for ctx just add all because it is source/partner already different entries
  merged_sv_bp_intra = get_range_as_bp_gr(merged_svs_gr[merged_svs_gr$svtype!="CTX"])
  
  merged_sv_bp_ctx = gr_midpoint(merged_svs_gr[merged_svs_gr$svtype=="CTX"])
  merged_sv_bp_ctx$sv_bp_orientation="ctx"
  
  merged_sv_bp_gr = c(merged_sv_bp_intra,merged_sv_bp_ctx)
  return(merged_sv_bp_gr)
}


get_recurrent_sv_overlaps = function(merged_svs_gr,svs_df,id_col="patient_sv_merged",svs_cols=NULL,svs_cols_display=NULL,flag_exclude_ctx=TRUE,
                                     reciprocal_overlap = 0.5, svtype_matching = FALSE,flag_within_tumor=FALSE) {
  
  if(length(merged_svs_gr[merged_svs_gr$svtype=="CTX"])>0 & flag_exclude_ctx){
    print("Warning, CTX svs present and will be ignored")
  } 
  if(flag_exclude_ctx){
    intra_sv_gr=merged_svs_gr[merged_svs_gr$svtype!="CTX"]
  } else {
    intra_sv_gr=merged_svs_gr
  }
  
  if(is.null(svs_cols)) {
    svs_cols = c(id_col,"patient_label","cytoband","tumor_af_mean","svlen_mean","cancer_type")
  }
  svs_cols = svs_cols[svs_cols %in% names(svs_df)]
  if(is.null(svs_cols_display)) {
    svs_cols_display = svs_cols[!svs_cols %in% c(id_col)]
  }
  
  overlaps = get_reciprocal_overlap_pairs(intra_sv_gr,intra_sv_gr,reciprocal_overlap = reciprocal_overlap, svtype_matching = svtype_matching,ignore_strand = T)
  overlaps = overlaps %>% 
    left_join(svs_df[,svs_cols],by=c("set1"=id_col)) %>% dplyr::rename_with(.cols=svs_cols_display,.fn=function(x){paste0("set1_",x)}) %>%
    left_join(svs_df[,svs_cols],by=c("set2"=id_col)) %>% dplyr::rename_with(.cols=svs_cols_display,.fn=function(x){paste0("set2_",x)}) 
  
  if(flag_within_tumor==T) {
    overlaps = overlaps %>% dplyr::filter(set1_patient_label==set2_patient_label & set1 != set2)
  } else {
    overlaps = overlaps %>% dplyr::filter(set1_patient_label!=set2_patient_label)
  }
  
  return(overlaps)
}


get_recurrent_ctx_overlaps = function(merged_svs_gr,svs_df,svs_cols=NULL,svs_cols_display=NULL,ctx_distance=1e6,
                                      reciprocal_overlap = 0.5, svtype_matching = FALSE,flag_within_tumor=FALSE) {
  if(length(merged_svs_gr[merged_svs_gr$svtype!="CTX"])>0){
    print("Warning, non CTX svs present and will be ignored")
  }
  
  if(is.null(svs_cols)) {
    svs_cols = c("patient_sv_merged","sv_merged_coordinate","start_gene_name","end_gene_name","patient_label","partner_sv_merged","partner_sv_merged_coordinate","cancer_type")
  }
  svs_cols = svs_cols[svs_cols %in% names(svs_df)]
  if(is.null(svs_cols_display)) {
    svs_cols_display = svs_cols[!svs_cols %in% c("patient_sv_merged")]
  }
  
  ctx_gr=merged_svs_gr[merged_svs_gr$svtype=="CTX"]
  
  ctx_overlaps = get_reciprocal_overlap_pairs(resize_gr_distance(ctx_gr,ctx_distance),resize_gr_distance(ctx_gr,ctx_distance),reciprocal_overlap = reciprocal_overlap, svtype_matching = svtype_matching,ignore_strand = T)
  
  ctx_overlaps = ctx_overlaps %>% 
    left_join(svs_df[,svs_cols],by=c("set1"="patient_sv_merged")) %>% dplyr::rename_with(.cols=svs_cols_display,.fn=function(x){paste0("set1_",x)}) %>%
    left_join(svs_df[,svs_cols],by=c("set2"="patient_sv_merged")) %>% dplyr::rename_with(.cols=svs_cols_display,.fn=function(x){paste0("set2_",x)}) 
  
  if(flag_within_tumor==T) {
    ctx_overlaps = ctx_overlaps %>% filter(set1_patient_label==set2_patient_label  & set1 != set2) %>% unique()
  }  else {   
    ctx_overlaps = ctx_overlaps %>% filter(set1_patient_label!=set2_patient_label) %>% unique()
  }
  #partner no match
  ctx_no_partner_match = ctx_overlaps %>% anti_join(ctx_overlaps[,c("set1_partner_sv_merged","set2_partner_sv_merged")],by=c("set1"="set1_partner_sv_merged","set2"="set2_partner_sv_merged"))
  #ctx_no_partner_match %>% filter(set1_start_gene_name=="EWSR1" & set2_end_gene_name=="FLI1") 
  
  #partner match
  ctx_overlaps = ctx_overlaps %>% inner_join(ctx_overlaps[,c("set1","set2","set1_partner_sv_merged","set2_partner_sv_merged")],
                                             by=c("set1"="set1_partner_sv_merged","set2"="set2_partner_sv_merged","set1_partner_sv_merged"="set1","set2_partner_sv_merged"="set2"))
  
  
  return(ctx_overlaps)
}


get_sv_bp_overlaps_ctx_self_other = function(merged_svs,sv_bp_cols,id_col="patient_sv_merged",partner_id_col="partner_sv_merged") {
  sv_bp_cols=unique(c(id_col,sv_bp_cols))
  sv_bp_overlaps_ctx_self = merged_svs %>%
    filter(svtype=="CTX") %>% select(all_of(c(sv_bp_cols,"svtype"))) %>%
    dplyr::rename_with(.cols=all_of(c(sv_bp_cols,"svtype")), .fn=function(x){paste0("set1_",x)}) %>%
    dplyr::rename(set1:=!!sym(paste0("set1_",id_col)))
  
  sv_bp_overlaps_ctx_other = merged_svs %>% 
    filter(svtype=="CTX") %>% select(all_of(c(sv_bp_cols,"svtype"))) %>%
    dplyr::rename_with(.cols=all_of(c(sv_bp_cols,"svtype")), .fn=function(x){paste0("set2_",x)}) %>%
    dplyr::rename(set2:=!!sym(paste0("set2_",id_col)))
  
  sv_bp_overlaps_ctx = sv_bp_overlaps_ctx_self %>%
    left_join(sv_bp_overlaps_ctx_other,by=c("set1"=paste0("set2_",partner_id_col)))
  
  return(sv_bp_overlaps_ctx)
}


lst_str = function(charstring){
  return(toString(sort(unique(charstring))))
}


get_clusters = function(clustered_merged_svs,merged_svs_gr) {
  clustered_merged_svs$cluster=as.character(clustered_merged_svs$cluster)
  
  clusters_df = clustered_merged_svs %>% group_by(cluster) %>%
    summarize(svs_cnt=length(unique(patient_sv_merged)),
              patient_cnt=length(unique(patient_label)),
              svs_lst=toString(sort(unique(patient_sv_merged))),
              patient_lst=toString(sort(unique(patient_label))),
    )
  clusters_df = clusters_df %>% arrange(-svs_cnt)
  
  
  clustered_merged_svs_gr = merged_svs_gr
  mcols(clustered_merged_svs_gr)  = mcols(clustered_merged_svs_gr) %>% as.data.frame() %>% left_join(clustered_merged_svs)
  clustered_merged_svs_gr = clustered_merged_svs_gr[!is.na(clustered_merged_svs_gr$cluster)]
  
  #cluster consensus region 
  ## doesnt work yet for ctx because reduce keeps as separate ranges since they are not overlapping
  intra_clustered_merged_svs_gr = clustered_merged_svs_gr[clustered_merged_svs_gr$svtype!="CTX"]
  
  merged_svs_clusters = GenomicRanges::reduce(GenomicRanges::split(intra_clustered_merged_svs_gr,intra_clustered_merged_svs_gr$cluster),ignore.strand=T)
  merged_svs_clusters = merged_svs_clusters %>% unlist()
  merged_svs_clusters$cluster=names(merged_svs_clusters)
  #names(merged_svs_clusters)=NULL
  merged_svs_clusters_coord = merged_svs_clusters %>% as.data.frame()
  
  merged_svs_clusters_coord = merged_svs_clusters_coord %>% get_gr_coordinate(attr_name = "cluster_consensus_coord")
  merged_svs_clusters_coord$cluster_width=merged_svs_clusters_coord$width
  
  
  ctx_clustered_merged_svs = clustered_merged_svs_gr[clustered_merged_svs_gr$svtype=="CTX"] %>% as.data.frame()
  ctx_clusters_coord = ctx_clustered_merged_svs %>% group_by(cluster) %>% summarize(seqnames=unique(seqnames),start=min(start),end=max(end))
  ctx_clusters_coord = ctx_clusters_coord %>% get_gr_coordinate(attr_name = "cluster_consensus_coord") %>% dplyr::mutate(cluster_width=end-start)
  #because partners are required will be separate entry 
  
  clusters_df$cluster=as.character(clusters_df$cluster)
  
  coord_cols = c("cluster","cluster_consensus_coord","cluster_width","seqnames","start","end")
  clusters_coord = rbind(merged_svs_clusters_coord[,coord_cols],ctx_clusters_coord[,coord_cols])
  clusters_df = clusters_df %>% left_join(clusters_coord)
  
  return(clusters_df)
}


get_gr_distance = function(coord_lst1,coord_lst2,flag_as_breakpoints=F,return_full_df=F) {
  #can return both range distance and bp distance 
  
  if(length(coord_lst1)!=length(coord_lst2)) {
    print("Need same length")
    print(paste0(length(coord_lst1)," < length list 1, length list 2> ",length(coord_lst2)))
    return(data.frame)
  }
  
  if(isS4(coord_lst1)) {
    gr1=coord_lst1
    coord_lst1 = gr1 %>% as.data.frame() %>% get_gr_coordinate() %>% select(coordinate) %>% flatten_chr()
  } else {
    gr1=GRanges(coord_lst1)
  }
  if(isS4(coord_lst2)) {
    gr2=coord_lst2
    coord_lst2 = gr2 %>% as.data.frame() %>% get_gr_coordinate() %>% select(coordinate) %>% flatten_chr()
  } else {
    gr2=GRanges(coord_lst2)
  }
  if(flag_as_breakpoints){
    #minimum of gr1 end- gr2 start and gr2 start - gr1 end should make it insensitive to whether gr1 or gr2 is earlier
    distance = mapply(FUN = min,abs(end(gr1)-start(gr2)),abs(start(gr1)-end(gr2)),
                      abs(start(gr1)-start(gr2)),abs(end(gr1)-end(gr2)))
    #also take into account overlapping/within svs so need to test all coordinates
    
  } else {
    #as ranges
    distance = GenomicRanges::distance(gr1,gr2,ignore.strand=T)
  }
  
  if(return_full_df) {
    overview=data.frame(coord1=coord_lst1)
    overview = cbind(overview,gr1 %>% as.data.frame() %>% select(c("seqnames","start","end")) %>% dplyr::rename_with(function(x){paste0("coord1_",x)}))
    
    overview$coord2=coord_lst2
    overview = cbind(overview,gr2 %>% as.data.frame() %>% select(c("seqnames","start","end")) %>% dplyr::rename_with(function(x){paste0("coord2_",x)}))
    
    overview$distance=distance
    overview$flag_as_breakpoints=flag_as_breakpoints
    
    #overview %>% head() %>% View()
    return(overview)
  } else {
    return(distance)
  }
}

get_sv_bp_overlaps_distance = function(merged_svs_gr,merged_svs,sv_bp_distance=5e6, sv_bp_cols=NULL,flag_within_tumor=NULL,id_col="patient_sv_merged",partner_id_col="partner_sv_merged") {
#general function to replace get_complex_sv_bp_overlaps
  if(is.null(flag_within_tumor)) {
    stop("Please specify flag_within_tumor to return either svs within distance of same patient or across patients, filter col patient_label")
  }
  if(is.null(sv_bp_cols)) {
    sv_bp_cols= c("patient_sv_merged","patient_label","cancer_type","sv_merged_coordinate","partner_sv_merged","partner_sv_merged_coordinate")
  }
  sv_bp_cols = c(id_col,partner_id_col,sv_bp_cols) %>% unique()
  
  merged_svs_gr_ids = unique(sort(mcols(merged_svs_gr)[,id_col] %>% as.data.frame() %>% flatten_chr()))
  merged_svs_ids = unique(sort(merged_svs[,id_col] ))
  
  if(!all(merged_svs_gr_ids==merged_svs_ids)) {
    print(paste0("merged_svs_gr and merged_svs not same ",id_col,", making intersection ..."))
    
    intersection = intersect(merged_svs_gr_ids,merged_svs_ids)
    print(paste0("gr: ",cnt_str(merged_svs_gr_ids), " df: ",cnt_str(merged_svs_ids), " intersection: ",cnt_str(intersection)))
    if(all(names(merged_svs_gr)==mcols(merged_svs_gr)[,id_col])) {
      merged_svs_gr = merged_svs_gr[intersection]
    } else {
      stop("Cannot subset merged_svs_gr, id col is not row name")
    }
    merged_svs = merged_svs %>% filter(!!sym(id_col) %in% intersection)
  }
  #intra svs split in start and end, for ctx just add all because it is source/partner already different entries
  merged_sv_bp_gr = get_svs_as_bp_gr(merged_svs_gr)
  
  sv_bp_overlaps = get_recurrent_sv_overlaps(resize_gr_distance(merged_sv_bp_gr,sv_bp_distance),svs_df=merged_svs,id_col=id_col,
                                             flag_exclude_ctx = F,flag_within_tumor = flag_within_tumor,
                                             reciprocal_overlap = 0,svtype_matching = F,svs_cols = sv_bp_cols)
  
  ## todo coordinate col as variable
  sv_bp_overlaps$distance=get_gr_distance(sv_bp_overlaps$set1_sv_merged_coordinate,sv_bp_overlaps$set2_sv_merged_coordinate)
  
  sv_bp_overlaps = sv_bp_overlaps %>% filter(distance<sv_bp_distance)
  
  #for all ctxes add their self/partner to the df 
  #they will never be in there normally so do not have to check
  #but remove the ones only consisting of partners
  ## that is done later now:
  ## remove clusters that are just edge ++ partner
  ##no_complex_sv = complex_sv_stats %>% filter(svs_lst==partner_svs_lst & svs_cnt==2)
  ##nrow(no_complex_sv)
  
  sv_bp_overlaps_ctx = get_sv_bp_overlaps_ctx_self_other(merged_svs,sv_bp_cols=sv_bp_cols,id_col=id_col,partner_id_col=partner_id_col)
  
  sv_bp_overlaps = rbind_no_colmatch(sv_bp_overlaps,sv_bp_overlaps_ctx)
  
  return(sv_bp_overlaps)
  
}

get_complex_sv_bp_overlaps = function(merged_svs_gr,merged_svs,complex_sv_bp_distance=5e6, complex_sv_bp_cols=NULL) {
  #complex_sv_bp_overlaps = get_sv_bp_overlaps_distance(merged_svs_gr,merged_svs,sv_bp_distance=complex_sv_bp_distance, sv_bp_cols=complex_sv_bp_cols,flag_within_tumor = T)
  #return(complex_sv_bp_overlaps)
  
  
  if(is.null(complex_sv_bp_cols)) {
    complex_sv_bp_cols= c("patient_sv_merged","patient_label","cancer_type","sv_merged_coordinate","partner_sv_merged","partner_sv_merged_coordinate")
  }
  #intra svs split in start and end, for ctx just add all because it is source/partner already different entries
  merged_sv_bp_gr = get_svs_as_bp_gr(merged_svs_gr)
  
  complex_sv_bp_overlaps = get_recurrent_sv_overlaps(resize_gr_distance(merged_sv_bp_gr,complex_sv_bp_distance),
                                                     svs_df=merged_svs,flag_exclude_ctx = F,flag_within_tumor = T,reciprocal_overlap = 0,svtype_matching = F,svs_cols = complex_sv_bp_cols)
  
  complex_sv_bp_overlaps$distance=get_gr_distance(complex_sv_bp_overlaps$set1_sv_merged_coordinate,complex_sv_bp_overlaps$set2_sv_merged_coordinate)
  
  complex_sv_bp_overlaps = complex_sv_bp_overlaps %>% filter(distance<complex_sv_bp_distance)
  
  #for all ctxes add their self/partner to the df 
  #they will never be in there normally so do not have to check
  #but remove the ones only consisting of partners
  ## that is done later now:
  ## remove clusters that are just edge ++ partner
  ##no_complex_sv = complex_sv_stats %>% filter(svs_lst==partner_svs_lst & svs_cnt==2)
  ##nrow(no_complex_sv)
  
  complex_sv_bp_overlaps_ctx = get_sv_bp_overlaps_ctx_self_other(merged_svs %>% filter(patient_sv_merged %in% merged_sv_bp_gr$patient_sv_merged),sv_bp_cols=complex_sv_bp_cols)
  
  complex_sv_bp_overlaps = rbind_no_colmatch(complex_sv_bp_overlaps,complex_sv_bp_overlaps_ctx)
  
  return(complex_sv_bp_overlaps)
 
}
get_complex_sv_bp_cluster_mapping = function(complex_sv_bp_overlaps,attr_name="complex_sv_cluster") {
  complex_sv_bp_graph = graph_from_data_frame(complex_sv_bp_overlaps[,c("set1","set2")],directed=F)
  complex_sv_bp_cluster_mapping = get_sv_cluster_mapping(complex_sv_bp_graph) %>% dplyr::rename(!!sym(attr_name):=cluster)
  return(complex_sv_bp_cluster_mapping)
}

get_ctx_incoming_outgoing = function(target_complex_sv_bp_overlaps,target_complex_sv_id=NULL) {
  #logical if set1/2 same chrom and ctx because these are physically nearby breakpoints
  #incoming and outgoing ctxes per chrom are mismatch between set1/set2
  # both ctx in and out are the chrom that occur in set1_chrom  && set2_chrom
  #target_complex_sv_bp_overlaps  %>% filter(set1_chrom!=set2_chrom)
  
  if(!is.null(target_complex_sv_id)) {
    target_complex_sv_bp_overlaps =  target_complex_sv_bp_overlaps %>% filter(complex_sv_id==target_complex_sv_id)
  }
  if(target_complex_sv_bp_overlaps$complex_sv_id %>% unique() %>% length() > 1) {
    print("Only safe for single complex, exiting")
    return(data.frame())
  }
  
  
  #should look at if all partner chrom are in the set1/set2
  #set1/set2 svs are "in" look at their partners "out" if these chrom are also in set 1/2
  #select 
  
  #select ctx and their partner chrom
  if(FALSE) {
    ctx_incoming_outgoing = merged_svs %>% filter(patient_sv_merged %in%
                                                    c(target_complex_sv_bp_overlaps$set1,target_complex_sv_bp_overlaps$set2)) %>% 
      filter(svtype=="CTX") %>% select(patient_sv_merged,chrom,partner_chrom)
  }
  
  #remove cases where set1_partner_sv_merged == set2
  target_complex_sv_bp_overlaps = target_complex_sv_bp_overlaps %>% 
    filter(set1_partner_sv_merged!=set2) 
  
  #either annotated or like above
  ctx_incoming_outgoing = target_complex_sv_bp_overlaps %>% 
    filter(set1_svtype=="CTX") %>%
    select(set1,set1_chrom,set1_partner_chrom) %>%
    dplyr::rename(patient_sv_merged=set1,chrom=set1_chrom,partner_chrom=set1_partner_chrom)
  
  ctx_incoming_outgoing = rbind(ctx_incoming_outgoing,
                                target_complex_sv_bp_overlaps %>% 
                                  filter(set2_svtype=="CTX") %>%
                                  select(set2,set2_chrom,set2_partner_chrom) %>%
                                  dplyr::rename(patient_sv_merged=set2,chrom=set2_chrom,partner_chrom=set2_partner_chrom)
  ) %>% unique()    
  
  #filter partner chrom should be in target complex sv bp overlaps
  ctx_incoming_outgoing = ctx_incoming_outgoing %>% 
    mutate(partner_in_complex = partner_chrom %in% ctx_incoming_outgoing$chrom)
  
  
  return(ctx_incoming_outgoing)
}

get_singular_edges = function(complex_sv_bp_overlaps,complex_sv_bp_cluster_mapping,vertex_names=F) {
  if(vertex_names) {
    complex_sv_bp_overlaps = complex_sv_bp_overlaps %>% 
      dplyr::mutate(set1=paste0(set1_chrom,"_",stri_sub(set1,-10,-1)),
                    set2=paste0(set2_chrom,"_",stri_sub(set2,-10,-1)))
  }
  
  complex_sv_edges = complex_sv_bp_overlaps[,c("set1","set2","complex_sv_id")]
  complex_sv_edges = complex_sv_edges %>% rowwise() %>% dplyr::mutate(uq = paste0(sort(c(set1,set2)),collapse = "_"))
  complex_sv_edges = complex_sv_edges[!duplicated(complex_sv_edges$uq),] %>% as.data.frame()
  return(complex_sv_edges)
}

## todo incoming/outgoing ctx 
#full list of overlaps to/from ctx only 
get_biconnected_components = function(graph) {
  biconn = igraph::biconnected_components(graph)
  biconn_lengths=lapply(biconn$components,length)
  
  cycle_vertices=data.frame()
  for(i in 1:length(biconn$components)) {
    vertex_lst = biconn$components[i] %>% unlist() %>% names()
    vertex_df = data.frame(vertex = vertex_lst)
    vertex_df$cycle_id=i
    vertex_df$cycle_size=biconn_lengths[[i]]
    cycle_vertices = rbind(cycle_vertices,vertex_df)
  }
  cycle_vertices = cycle_vertices %>% dplyr::mutate(is_cycle=cycle_size>2)
  return(cycle_vertices)
}

annotate_complex_sv_bp_overlaps = function(complex_sv_bp_overlaps,complex_sv_bp_cluster_mapping,merged_svs) {
  complex_sv_bp_cluster_mapping[duplicated(complex_sv_bp_cluster_mapping$patient_sv_merged),]  %>% nrow() ==0
  
  complex_sv_bp_overlaps = complex_sv_bp_overlaps %>% #[,c("set1","set1_svtype","set2")] %>% 
    left_join(complex_sv_bp_cluster_mapping[,c("patient_sv_merged","complex_sv_id")],by=c("set1"="patient_sv_merged"),relationship = "many-to-many") %>%
    left_join(complex_sv_bp_cluster_mapping[,c("patient_sv_merged","complex_sv_id")],by=c("set2"="patient_sv_merged"),relationship = "many-to-many")
  
  if(complex_sv_bp_overlaps %>% filter(complex_sv_id.x!=complex_sv_id.y) %>% nrow() == 0) {
    complex_sv_bp_overlaps = complex_sv_bp_overlaps %>% select(-complex_sv_id.y) %>% dplyr::rename(complex_sv_id=complex_sv_id.x)
  }
  
  #necessary for for inferring cycles automated and plots
  complex_sv_bp_overlaps = complex_sv_bp_overlaps %>% left_join(merged_svs[,c("patient_sv_merged","chrom","partner_chrom","partner_sv_merged")],by=c("set1"="patient_sv_merged")) %>% dplyr::rename(set1_chrom=chrom,set1_partner_chrom=partner_chrom) %>% #set1_partner_sv_merged=partner_sv_merged) %>%
    left_join(merged_svs[,c("patient_sv_merged","svtype","chrom","partner_chrom")],by=c("set2"="patient_sv_merged")) %>% dplyr::rename(set2_chrom=chrom,set2_partner_chrom=partner_chrom)
  
  return(complex_sv_bp_overlaps)
}


is_complex_sv_cyclical = function(complex_sv_bp_overlaps,target_complex_sv_id,return_cycles=F,plot=F,verbose=F,strict=T) {
  #target_complex_sv_id = "M552AAA_complex_ff3072acf71410ff81975bf9854bc9cf"
  target_complex_sv_bp_overlaps = complex_sv_bp_overlaps %>% filter(complex_sv_id==target_complex_sv_id)
  
  if(!"set1_chrom" %in% colnames(target_complex_sv_bp_overlaps)) {
    print("Could not check incoming/outgoing ctxes, chrom info missing")
  } else {
    ctx_incoming_outgoing = get_ctx_incoming_outgoing(target_complex_sv_bp_overlaps)
    
    #are there ctxes with chromosomes not in the clustered svs?
    ctx_loose_ends = ctx_incoming_outgoing %>% filter(partner_in_complex==F)
    
    if(ctx_loose_ends %>% nrow() > 0) {
      if(strict) { return(F) } 
      if(verbose) { 
        print(paste0(target_complex_sv_id, ": chromosomes present without both incoming and outgoing ctxes"))  
        print("Filtering out")
        print(ctx_loose_ends) 
      }
      target_complex_sv_bp_overlaps = target_complex_sv_bp_overlaps %>% 
        filter(! set1 %in% ctx_loose_ends$patient_sv_merged & !set2 %in% ctx_loose_ends$patient_sv_merged)
      
    } else if(verbose) {
      print(paste0(target_complex_sv_id, ": ctx_incoming_outgoing pass"))
      print(ctx_incoming_outgoing)
    }
  }
  
  
  target_complex_sv_edges = get_singular_edges(target_complex_sv_bp_overlaps)
  
  if(target_complex_sv_edges %>% nrow()==0) {return(F)}
  
  if(plot) {
    target_complex_sv_edges_plot = get_singular_edges(target_complex_sv_bp_overlaps,vertex_names = T)
    target_complex_sv_bp_graph = graph_from_data_frame(target_complex_sv_edges_plot[,c("set1","set2")],directed=F)
    #if(plot) {  plot(target_complex_sv_bp_graph,vertex.labels=NA) %>% print() }
    plot(target_complex_sv_bp_graph,main=target_complex_sv_id,size=5)  %>% print() 
  }
  
  target_complex_sv_bp_graph = graph_from_data_frame(target_complex_sv_edges[,c("set1","set2")],directed=F)
  
  
  #igraph::eulerian_cycle(target_complex_sv_bp_graph) ## too strict 
  #igraph::eulerian_path(target_complex_sv_bp_graph)
  
  cycle_vertices = get_biconnected_components(target_complex_sv_bp_graph)
  vertices_not_in_cycle = target_complex_sv_edges %>% filter(!set1 %in% filter(cycle_vertices,is_cycle)$vertex)
  if(vertices_not_in_cycle %>% nrow() > 0) {
    vertices_not_in_cycle$vertex = vertices_not_in_cycle$set1
    vertices_not_in_cycle$is_cycle = F
    cycle_vertices=rbind_no_colmatch(cycle_vertices,vertices_not_in_cycle[,c("vertex","is_cycle")])
  }
  
  if(return_cycles) {
    return(cycle_vertices)
  } else {
    if(strict) {
      #true if perfectly connected/true cyclical
      # todo consider requiring single cycle
      return((vertices_not_in_cycle %>% nrow() == 0))
    } else {
      return((filter(cycle_vertices,cycle_size>3) %>% nrow() > 0))
    }
  }
}
annotate_complex_sv_stats_cycles = function(complex_sv_stats,complex_sv_bp_overlaps,get_cycles_global=T,plot=T) {
  complex_sv_stats$flag_cycle=F
  complex_sv_stats$flag_cycle_present=F
  complex_sv_stats$flag_cycle_ctx_only=F
  complex_sv_stats$flag_cycle_cnt=0
  if(get_cycles_global) { cohort_cycles=data.frame() }
  
  
  for(target in complex_sv_stats$complex_sv_id) {
    #at least one ctx otherwise just looking at overlapping intra svs
    
    target_ctx_overlaps = complex_sv_bp_overlaps %>% filter(complex_sv_id == target) %>% 
      filter(set1_svtype=="CTX" | set2_svtype=="CTX")
    
    if(nrow(target_ctx_overlaps) == 0) { next() }
    complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle_present")] = is_complex_sv_cyclical(complex_sv_bp_overlaps,target,verbose=F,strict=F)
    
    complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle")] = is_complex_sv_cyclical(complex_sv_bp_overlaps,target,verbose=F)
    
    target_ctx_only_overlaps = target_ctx_overlaps %>% 
      filter(set1_svtype=="CTX" & set2_svtype=="CTX")
    
    complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle_ctx_only")] = is_complex_sv_cyclical(target_ctx_only_overlaps,target,verbose=F)
    
    if(complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle")] == T) {
      cycles =  is_complex_sv_cyclical(complex_sv_bp_overlaps,target,verbose=F,plot = T,return_cycles = T)
      cycle_cnt=length(unique(cycles$cycle_id))
      cycles$complex_sv_id=target
      cycles$strict=T
      cycles$cycle_cnt=cycle_cnt
      if(get_cycles_global) { cohort_cycles=rbind(cohort_cycles,cycles) }
      
      complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle_cnt")] = unique(cycles$cycle_cnt)
      
      if(plot){
        pdf(paste0(cohort_results_dir,"complex_sv_clusters/complex_cycles.",target,".pdf"))
        is_complex_sv_cyclical(complex_sv_bp_overlaps,target,verbose=F,plot = T)
        dev.off()
      }
      
    } else if(complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle_present")] == T) {
      #strict == false
      # todo refactor this
      cycles =  is_complex_sv_cyclical(complex_sv_bp_overlaps,target,verbose=F,plot = T,return_cycles = T,strict=F)
      cycle_cnt=length(unique(cycles$cycle_id))
      
      cycles$complex_sv_id=target
      cycles$strict=F
      cycles$cycle_cnt=cycle_cnt
      if(get_cycles_global) { cohort_cycles=rbind(cohort_cycles,cycles) }
      
      complex_sv_stats[complex_sv_stats$complex_sv_id==target,c("flag_cycle_cnt")] = cycle_cnt
      
      
      if(plot){
        pdf(paste0(cohort_results_dir,"complex_sv_clusters/complex_cycles.partial.",target,".pdf"))
        is_complex_sv_cyclical(complex_sv_bp_overlaps,target,verbose=F,plot = T,strict=F)
        dev.off()
      }
    }
    
    if(get_cycles_global) { cohort_cycles <<- cohort_cycles }
    
  }
  return(complex_sv_stats)
  
}


### Peak calling 
##also for cov svs so in general

call_peaks = function(reduced_gr,max_peaks_only=T) {
  cov_obj = reduced_gr %>% coverage()
  
  peak_collection=GRanges()
  for(chrom in names(cov_obj)[!grepl("_",names(cov_obj))]) {
    if(length(cov_obj[[chrom]])==0) {next()}
    
    max_peak= max(cov_obj[[chrom]])
    
    if(max_peaks_only==T) {
      peak_heights=c(max_peak, (max_peak-1))
    } else {
      peak_heights=1:max_peak
    }
    
    peak_collection_chrom=GRanges()
    
    for(peak_height in peak_heights) {
      cov_shared_region = IRanges::slice(cov_obj[[chrom]],lower=peak_height,upper = peak_height)
      if(length(cov_shared_region)==0) next()
      peaks = GRanges(IRanges(cov_shared_region),seqnames=chrom)
      peaks$cov = peak_height
      peak_collection_chrom=c(peak_collection_chrom,peaks)
    }
    
    peak_collection_chrom$peak_id = paste0("peak_",chrom,"_",1:length(peak_collection_chrom))
    names(peak_collection_chrom) = peak_collection_chrom$peak_id
    
    peak_collection=c(peak_collection,peak_collection_chrom)
  }
  
  peak_collection_df = as.data.frame(peak_collection)
  peak_collection_df = peak_collection_df %>% dplyr::rename(peak_width = width)
  peak_collection_df$peak_coordinate= paste0(peak_collection_df$seqnames,":",peak_collection_df$start,"-",peak_collection_df$end)
  
  return(peak_collection_df)
}

get_sv_peaks_overlaps = function(svs_gr,svs_df,peaks_df,svs_id_col="patient_sv_merged",peak_id_col="peak_id",
                                 sv_peak_cols=NULL,peak_cols=NULL){
  if(is.null(sv_peak_cols)) {
    sv_peak_cols = c("patient_sv_merged","patient_label","cytoband","tumor_af_mean","svlen_mean","cancer_type")
  }
  if(is.null(peak_cols)) {
    peak_cols=c("peak_id","seqnames","start","end","peak_coordinate")
  }
  
  peaks_gr=GRanges(peaks_df)
  names(peaks_gr)=peaks_df[,peak_id_col]
  
  sv_peaks_overlaps = get_reciprocal_overlap_pairs(svs_gr,peaks_gr,reciprocal_overlap = 0,svtype_matching = F)
  sv_peaks_overlaps = sv_peaks_overlaps %>% dplyr::rename(!!svs_id_col:=set1,!!peak_id_col:=set2) %>% 
    left_join(svs_df[,c(svs_id_col,sv_peak_cols)],by=svs_id_col) %>% left_join(peaks_df[,c(peak_id_col,peak_cols)],by=peak_id_col) 
  
  return(sv_peaks_overlaps) 
}

get_patients_per_peak = function(sv_peaks_overlaps,group_cols=c("peak_id"),no_cancer_type=F) {
  patients_per_peak = sv_peaks_overlaps  %>% 
    group_by(across(all_of(group_cols))) %>%
    summarize(patient_lst=lst_str(patient_label),
              patient_cnt=cnt_str(patient_label),
              patient_md5=md5(patient_lst),
              cancer_type_lst=lst_str(cancer_type),
              cancer_type_cnt=cnt_str(cancer_type)
    ) %>% arrange(-patient_cnt,-cancer_type_cnt)
  if(no_cancer_type) {
    patients_per_peak = sv_peaks_overlaps  %>% 
      group_by(across(all_of(group_cols))) %>%
      summarize(patient_lst=lst_str(patient_label),
                patient_cnt=cnt_str(patient_label),
                patient_md5=md5(patient_lst)
      ) %>% arrange(-patient_cnt)
  }
  return(patients_per_peak)
}

get_cancer_type_peaks = function(input_gr,peak_type_value,return_map=F) {
  cancer_type_sv_peaks = data.frame()
  cancer_type_map_sv_peaks_patients = data.frame()
  
  for(cancer_type_value in unique(cohort$cancer_type)) {
    patient_subset = filter(cohort,cancer_type==cancer_type_value)$patient_id
    
    input_gr_subset = input_gr[input_gr$patient_id %in% patient_subset]
    reduced_input_gr_subset  = unlist(GenomicRanges::reduce(split(input_gr_subset, input_gr_subset$patient_id),ignore.strand=T))
    reduced_input_gr_subset$patient_id=names(reduced_input_gr_subset)
    
    input_gr_subset_peaks = call_peaks(reduced_input_gr_subset,max_peaks_only = F)
    input_gr_subset_peaks$peak_type=peak_type_value
    
    input_gr_subset_peaks$cancer_type=cancer_type_value
    input_gr_subset_peaks$peak_id=paste0(input_gr_subset_peaks$cancer_type,"_",input_gr_subset_peaks$peak_id)
    
    #below not necessary
    #cancer_type_sv_peaks=rbind(cancer_type_sv_peaks,input_gr_subset_peaks)
    
    sv_peaks_overlaps = get_sv_peaks_overlaps(svs_gr = input_gr_subset,svs_df = cohort_svs_df,peaks_df = input_gr_subset_peaks,svs_id_col="patient_sv_merged",peak_id_col="peak_id") 
    sv_peaks_overlaps$cancer_type=cancer_type_value
    
    patients_per_peak = get_patients_per_peak(sv_peaks_overlaps,group_cols=c("peak_id"),no_cancer_type=T) 
    
    input_gr_subset_peaks = input_gr_subset_peaks %>% left_join(patients_per_peak)
    cancer_type_sv_peaks=rbind(cancer_type_sv_peaks,input_gr_subset_peaks)
    
    cancer_type_map_sv_peaks_patients = rbind(cancer_type_map_sv_peaks_patients, sv_peaks_overlaps %>% select(all_of(c(sv_peak_cols,peak_cols,"cancer_type"))) )
  }
  
  if(return_map) {
    return(cancer_type_map_sv_peaks_patients)
  }
  return(cancer_type_sv_peaks)
  
}

get_cn_peaks = function(segments,segments_df,cn_call_lst=c("gain","loss","loh"),return_map=F) {
  if(!exists("seg_cols")) {
    seg_cols=c("patient_cna_id","patient_label","call","cancer_type")
  }
  if(!exists("peak_cols")) {
    peak_cols=c("peak_id","seqnames","start","end","peak_coordinate")
  }
  cn_peaks = data.frame()
  map_cn_peaks_patients = data.frame()
  
  for(cn_type in cn_call_lst){
    #cn_type="gain"
    cohort_calls = segments[segments$call==cn_type]
    
    #remove blacklist
    #  overlap_blacklist = get_reciprocal_overlap_pairs(cohort_calls,blacklist,reciprocal_overlap = 0,svtype_matching = F)
    #  cohort_calls = cohort_calls[!names(cohort_calls) %in% overlap_blacklist$set1] 
    
    reduced_cohort_calls = unlist(GenomicRanges::reduce(split(cohort_calls, cohort_calls$patient_label),ignore.strand=T))
    reduced_cohort_calls$patient_label = names(reduced_cohort_calls)
    
    cn_coverage_df = call_peaks(reduced_cohort_calls,max_peaks_only = F)
    cn_coverage_df$cn_type=cn_type
    cn_coverage_df$peak_id=paste0(cn_coverage_df$peak_id,"_",cn_coverage_df$cn_type)
    
    cn_coverage_gr=GRanges(cn_coverage_df)
    names(cn_coverage_gr)=cn_coverage_gr$peak_id
    
    cn_peaks_overlaps = get_reciprocal_overlap_pairs(cohort_calls,cn_coverage_gr,reciprocal_overlap = 0,svtype_matching = F)
    cn_peaks_overlaps = cn_peaks_overlaps %>% dplyr::rename(patient_cna_id=set1,peak_id=set2)
    
    cn_peaks_overlaps = cn_peaks_overlaps %>% left_join(segments_df[,seg_cols],by="patient_cna_id") %>% left_join(cn_coverage_df[,peak_cols],by="peak_id") 
    
    map_cn_peaks_patients = rbind(map_cn_peaks_patients, cn_peaks_overlaps%>% select(all_of(c(seg_cols,peak_cols))) )
    
    patients_per_peak = get_patients_per_peak(cn_peaks_overlaps)
    # patients_per_peak = cn_peaks_overlaps  %>% 
    #   group_by(peak_id) %>% 
    #   summarize(patient_lst=toString(unique(sort(patient_id))),
    #             patient_cnt=length(unique(patient_id)),
    #             patient_md5=md5(patient_lst))
    
    cn_coverage_df = cn_coverage_df %>% left_join(patients_per_peak)
    cn_peaks=rbind(cn_peaks,cn_coverage_df)
    
    
    ## bins, deprecated
    # cn_cov_by_bins_entry = get_coverage_by_bins(cn_coverage_df,genome_bins)
    # cn_cov_by_bins_entry$peak_types=paste0("cn_",cn_type)
    # cn_cov_by_bins = rbind(cn_cov_by_bins,cn_cov_by_bins_entry)
  }
  if(return_map==T) {
    return(map_cn_peaks_patients)
  }
  return(cn_peaks)
}



## SV and SV bp distance functions ----
## from complex svs rmd

get_pairwise_distances = function(gr1, gr2 = NULL, ignore.strand = FALSE, ...) {
  library(reshape2)
  
  if (is.null(gr2)){
    gr2 = gr1
  }
  
  if (ignore.strand){
    strand(gr1) = '*'
    strand(gr2) = '*'
  }
  
  ix1 = rep(1:length(gr1), length(gr2))
  ix2 = rep(1:length(gr2), each = length(gr1))
  
  out = matrix(suppressWarnings(distance(gr1[ix1], gr2[ix2], ...)), nrow = length(gr1), ncol = length(gr2))
  out = out %>% melt()
  colnames(out) = c("lead_number","other","distance")
  return(out)
}

#set = merged_sv_bp_patients_grl["M002AAB"][[1]]
#set=merged_sv_bp_complex_grl["M226AAD_complex_8ace4e6f7d18e2059c91ab793b6aeb3f"][[1]]

get_chrom_lead_numbers = function(set,attr="lead_number") {
  if(isS4(set)) { 
    set=set[order(set@seqnames,set@ranges@start)]
    mcols(set)[,attr] = 1:length(set)
  } else {
    set=set %>% arrange(seqnames,start)
    set[,attr] = 1:nrow(set)
  }
  return(set)
}

annotate_svs_distance_to = function(set,id_col="patient_sv_merged",analyse_sv_bp=T,analysis_type="nearest") {
  set=set %>% get_chrom_lead_numbers()
  set_anno = mcols(set) %>% as.data.frame() 
  rownames(set_anno)=NULL
  pairwise_distances =  get_pairwise_distances(set,ignore.strand = T)
  pairwise_distances = pairwise_distances %>% filter(!is.na(distance) & lead_number!=other)
  pairwise_distances = pairwise_distances %>% left_join(set_anno[,c("lead_number",id_col)],by="lead_number") %>% left_join(set_anno[,c("lead_number",id_col)],by=c("other"="lead_number")) 
  
  pairwise_distances = pairwise_distances %>% dplyr::rename(!! sym(id_col)  := !! sym(paste0(id_col,".x")))
  
  #selected nearest, NA is dropped as well
  if(analysis_type=="nearest") {
    #make sure not self if bps
    
    pairwise_distances =  pairwise_distances %>% filter( !! sym(id_col) != !! sym(paste0(id_col,".y")))
    
    #note 2023-04-20 this is the original implementation, not sure if needed
    if(analyse_sv_bp) {  
      minimize_col=id_col
      #using id col -> even if started with bps, reporting smallest per sv
      #so then the smallest dist for either bp is chosen
    } else {
      #ranges
      minimize_col="lead_number" #otherwise put id_col
    }
    
    pairwise_distances = pairwise_distances[pairwise_distances$distance == ave(pairwise_distances$distance, pairwise_distances[,minimize_col],FUN=min),]
    pairwise_distances =   pairwise_distances %>% dplyr::rename(distance_to_nearest=distance) %>% 
      group_by(!!sym(minimize_col),distance_to_nearest) %>% 
      summarize(!! sym(paste0("nearest_",id_col)) := toString( unique(sort( !! sym(paste0(id_col,".y"))))))
  } 
  
  #consecutive remove overlapping that have distance == 0 then select nearest from the ones with higher/lower lead number
  #NB: all values in upstream col will also occur in downstream col so you can select either for postprocessing
  #for ranges: overlapping are removed 
  #for breakpoints: as dsDNA breaks regardless if same SV or not 
  #use ranges eg for amplicon CN segments and breakpoints for svs
  if(analysis_type=="consecutive") {
    #make sure to do ave on lead number to keep results for each breakpoint
    minimize_col="lead_number" #otherwise put id_col
    
    distances_downstream = pairwise_distances %>% filter(other>lead_number & distance>0)
    distances_downstream = distances_downstream[distances_downstream$distance == ave(distances_downstream$distance, distances_downstream[,minimize_col],FUN=min),]
    
    distances_downstream = distances_downstream %>% dplyr::rename(downstream_distance=distance) %>%
      group_by(!!sym(minimize_col),downstream_distance) %>%
      summarize(downstream_lead_number=lst_str(other),
                !! sym(paste0("downstream_",id_col)) := lst_str(!! sym(paste0(id_col,".y"))))
    
    distances_upstream = pairwise_distances %>% filter(other<lead_number & distance>0)
    distances_upstream = distances_upstream[distances_upstream$distance == ave(distances_upstream$distance, distances_upstream[,minimize_col],FUN=min),]
    
    distances_upstream = distances_upstream %>% dplyr::rename(upstream_distance=distance) %>%
      group_by(!!sym(minimize_col),upstream_distance) %>%
      summarize(upstream_lead_number=lst_str(other),
                !! sym(paste0("upstream_",id_col)) := lst_str(!! sym(paste0(id_col,".y"))))
    
    distances_updownstream = pairwise_distances[,c(id_col,"lead_number")] %>% unique() %>% left_join(distances_upstream)  %>% left_join(distances_downstream)
    pairwise_distances = distances_updownstream 
    
  }
  
  set_anno = set_anno %>% left_join(pairwise_distances)
  
  mcols(set) = set_anno
  return(set)
}

get_distance_to_df = function(split_grl,distance_to_nearest_id_col=c("patient_sv_merged"),id_cols=c("patient_sv_merged","complex_sv_id"),analyse_sv_bp=T,analysis_type="nearest") {
  
  split_grl = lapply(split_grl,function(x) { annotate_svs_distance_to(set=x,id_col = distance_to_nearest_id_col,analyse_sv_bp = analyse_sv_bp,analysis_type=analysis_type)})
  split_gr= split_grl[[names(split_grl)[1]]]
  for(i in names(split_grl)){
    split_gr = c(split_gr,split_grl[[i]])
  }
  distance_to_nearest_df = mcols(split_gr) %>% as.data.frame() #%>% select(all_of(id_cols),lead_number,paste0("nearest_",distance_to_nearest_id_col),distance_to_nearest) 
  rownames(distance_to_nearest_df) = NULL
  
  return(distance_to_nearest_df)
}


## CN/SV integration ----
get_chrom_cn_fractions = function(segments,segments_df,ranges,ranges_df,cohort,ranges_id_col="chr_arm",ranges_width_col="chr_arm_width",
                                  relative_cr=T,return_wide=T,add_chrom_summary=T,cna_id_col="patient_cna_id") {
  #can use it as this for chromosomes
  # get_chrom_cn_fractions(segments,segments_df,chromosomes,chromosomes_df,cohort,ranges_id_col="seqnames",ranges_width_col="chrom_width")
  
  #baseline shift already in segments
  chrom_cna_overlaps = get_reciprocal_overlap_pairs(segments,ranges,reciprocal_overlap = 0,svtype_matching = F)
  chrom_cna_overlaps = chrom_cna_overlaps %>% dplyr::rename(!!cna_id_col:=set1,!!ranges_id_col:=set2)
  
  modeled_seg_cols=c(cna_id_col,"width","cr_l2fc_50","maf_50")
  chrom_cna_overlaps = chrom_cna_overlaps %>% left_join(segments_df[,c("patient_label",modeled_seg_cols)] %>% unique(),by=c(cna_id_col))
  chrom_cna_overlaps = chrom_cna_overlaps %>% call_cna(cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4)
  
  #make relative copy ratio
  if(relative_cr) {
    chrom_cna_overlaps = chrom_cna_overlaps  %>% left_join(cohort[,c("patient_label","tumor_id")]) %>% get_relative_cr_l2fc(cohort) %>% rename_relative_cr_to_normal()
  }
  chrom_cn_summary_cols = c("frac_covered","cr_l2fc_50","relative_seg_covered_cr_stable","seg_covered")
  if("patient_label" %in% names(segments_df) & "patient_label" %in% names(ranges_df)) {
    #merg ensures filtering by match of patient label
    chrom_cna_overlaps = chrom_cna_overlaps %>% merge(ranges_df[,c(ranges_id_col,ranges_width_col,"patient_label")],by=c("patient_label",ranges_id_col)) %>% unique()
  } else {
    chrom_cna_overlaps = chrom_cna_overlaps %>% left_join(ranges_df[,c(ranges_id_col,ranges_width_col)],by=c(ranges_id_col)) %>% unique()
  }
  chrom_cn_fractions = call_cn_range_fraction(chrom_cna_overlaps,range_length_col=ranges_width_col,group_cols=c("patient_label",ranges_id_col))
  
  #need to divide by fraction to know whether all gain segments are within certain value relative to all gains instead of entire window.
  chrom_cn_fractions = call_cn_seg_same_state(chrom_cna_overlaps,chrom_cn_fractions,group_cols=c("patient_label",ranges_id_col,"call")) %>%
    dplyr::mutate(relative_seg_covered_cr_stable=seg_covered_cr_stable/frac_covered)
  if(return_wide==F) { return(chrom_cn_fractions) }
  
  chrom_cn_wide = chrom_cn_fractions %>% 
    dplyr::select(all_of(c("patient_label",ranges_id_col,"call",chrom_cn_summary_cols))) %>% 
    pivot_wider(names_from="call",values_from=all_of(chrom_cn_summary_cols) )
  
  if(return_wide==T & add_chrom_summary==T) {
    cohort_chrom_cna = get_overlaps_modeled_seg_summary_short(chrom_cna_overlaps,group_cols=c("patient_label",ranges_id_col)) %>% call_cna() 
    cohort_chrom_cna = call_cn_seg_same_state(chrom_cna_overlaps,cohort_chrom_cna,group_cols=c("patient_label",ranges_id_col))
    chrom_cn_wide = chrom_cn_wide %>% left_join(cohort_chrom_cna[,c("patient_label",ranges_id_col,"call","seg_covered_cr_data","cr_l2fc_50","maf_50","seg_covered_cr_stable")])
  }
  
  return(chrom_cn_wide)
}

get_chrom_centered_cn_segments = function(segments,segments_df,chromosomes,chromosomes_df,cohort) {
  
  chrom_cn_fractions = get_chrom_cn_fractions(segments,segments_df,
                                              chromosomes,chromosomes_df,cohort,ranges_id_col="seqnames",ranges_width_col="chrom_width",
                                              relative_cr = F,return_wide = T)
  chrom_cn_fractions = chrom_cn_fractions %>% dplyr::mutate(chrom_mean_cr=cr_l2fc_50)  
  chromosome_cn_calls = chrom_cn_fractions %>% select(patient_label,seqnames,chrom_mean_cr) %>% ungroup()
  
  chrom_centered_segments_df = segments_df %>% left_join(chromosome_cn_calls) %>% 
    dplyr::mutate(cr_l2fc_50_uncorrected=cr_l2fc_50,cr_l2fc_50=(cr_l2fc_50-chrom_mean_cr))
  
  chrom_centered_segments_df = chrom_centered_segments_df %>% call_cna()
  chrom_centered_segments_df = chrom_centered_segments_df %>% get_gr_coordinate("cna_coordinate")
  
  
  return(chrom_centered_segments_df)
  
}

get_cna_merged_gainloss_per_patient = function(split_segments,segments_df,cnas_id_col = "patient_cna_id",group_cols="",
                                               cnas_df_cols =  NULL,gap_size=6e3,split_by_cntype=c("gain","loss")) {
  #aka get_amplicons_merged_seg()
  #gap_size=6e3 default 
  if(is.null(cnas_df_cols)) {
    cnas_df_cols = c("patient_cna_id","cr_l2fc_50","maf_50","call","cna_coordinate")
  }
  cnas_df_cols=c(cnas_df_cols,group_cols) %>% unique()
  
  split_segments_grl = GenomicRanges::split(split_segments,split_segments$patient_label)
  split_segments_grl = lapply(split_segments_grl,function(x) { 
    get_cna_merged_gainloss(x, segments_df,cnas_id_col = cnas_id_col,cnas_df_cols =  cnas_df_cols,
                            group_cols=group_cols,gap_size=gap_size,split_by_cntype=split_by_cntype)})
  
   
  merged_seg = GRanges()
  for(i in names(split_segments_grl)){
    patient_split_segments_gr = split_segments_grl[[i]]
    patient_split_segments_gr$patient_label=i
    names(patient_split_segments_gr) = paste0(i,"_",names(patient_split_segments_gr))
    merged_seg = c(merged_seg,patient_split_segments_gr)
  }
  return(merged_seg)
}

#wrapper 
get_mega_merged_segments = function(segments,segments_df,gap_size=5e6,cnas_df_cols=NULL,group_cols=NULL,return_as_df=T,split_by_cntype=c("gain","loss")) {
  #main id col: patient_cna_id
  if(is.null(cnas_df_cols)) {
    cnas_df_cols = c("patient_cna_id","cr_l2fc_50","maf_50","call","patient_label")
  }
  if(is.null(group_cols)) {
    group_cols = c("patient_label")
  }
  
  split_by_cntype = split_by_cntype[split_by_cntype %in% unique(segments$call)]
  
  mega_merge_collection = GRanges()
  for(cntype in split_by_cntype)  {
    subset_chrom_segments = segments[segments$call==cntype]
  
    subset_chrom_segments_merged = get_cna_merged_gainloss_per_patient(
      subset_chrom_segments,segments_df,cnas_df_cols = cnas_df_cols,
      group_cols = group_cols,gap_size = gap_size,split_by_cntype = c(cntype))
    subset_chrom_segments_merged$patient_cna_id=names(subset_chrom_segments_merged)
    #needs to use whatever was used to split to prevent clashes in names
    subset_chrom_segments_merged$patient_cna_id = paste0(subset_chrom_segments_merged$patient_cna_id,"_",cntype)
    names(subset_chrom_segments_merged)=subset_chrom_segments_merged$patient_cna_id
    
    mega_merge_collection = c(mega_merge_collection,subset_chrom_segments_merged)
  }
  
  #can be mismatch between selection on call segments gr and df data intentionally
  #causes confusing cna merged names 
  mega_merge_collection = mega_merge_collection %>% call_cna()
  
  
  if(return_as_df){
    return(mega_merge_collection %>% as.data.frame())
  } else {
    return(mega_merge_collection)
  }
}


make_chromosome_cn_calls = function(segments,segments_df,chromosomes,chromosomes_df,chr_arms,chr_arms_df,make_relative_cr=F,
                                    threshold_chrom_frac_cn_state=0.7,threshold_chr_arm_frac_cn_state=0.7) {
  #first look at chromosome as a whole, then at chr arms
  
  #if subset of segments are used, still fractions relative to entire chromosome
  #so the fractions covered are NOT divided by the seg covered cr data 
  #example if seg covered cr data is 0.83 -> 
  #> 0.74/0.997 => 0.7422267
  #> 0.74/0.83 => 0.8915663 actually of the numerical part more is loh 
  #but easier to keep it like this so I can just select >0.7 as criterium and choose highest. 
  #the relative_seg_covered_cr_stable_loh is also NOT scaled so might be underestimation (but that is OK)
  
  chrom_cn_fractions = get_chrom_cn_fractions(segments,segments_df,
                                              chromosomes,chromosomes_df,cohort,ranges_id_col="seqnames",ranges_width_col="chrom_width",
                                              relative_cr = make_relative_cr,return_wide = F)
  
  
  chrom_cn_fractions$patient_chrom = paste0(chrom_cn_fractions$patient_label,"_",chrom_cn_fractions$seqnames)
  
  #select largest frac
  chrom_cn_fractions_highest = chrom_cn_fractions[
    chrom_cn_fractions$frac_covered == ave(chrom_cn_fractions$frac_covered, chrom_cn_fractions$patient_chrom,FUN=max),]
  chrom_cn_fractions_highest = chrom_cn_fractions_highest %>% dplyr::mutate(chrom_frac_stable = seg_covered_cr_stable > threshold_chrom_frac_cn_state) 
  chrom_cn_fractions_highest$region="chrom"
  chrom_cn_fractions_highest = chrom_cn_fractions_highest %>% dplyr::mutate(selected=chrom_frac_stable) %>% as.data.frame()
  
  
  ## chr arm
  
  chr_arm_cn_fractions = get_chrom_cn_fractions(segments,segments_df,chr_arms,chr_arms_df,cohort,relative_cr = make_relative_cr,return_wide = F)
  #chr_arm_cn_fractions = chr_arm_cn_fractions %>% filter(!grepl("acen|gvar|stalk",chr_arm))
  chr_arm_cn_fractions = chr_arm_cn_fractions %>% left_join(chr_arms_df[,c("chr_arm","seqnames")]) 
  
  chr_arm_cn_fractions$patient_chr_arm = paste0(chr_arm_cn_fractions$patient_label,"_",chr_arm_cn_fractions$chr_arm)
  chr_arm_cn_fractions$patient_chrom = paste0(chr_arm_cn_fractions$patient_label,"_",chr_arm_cn_fractions$seqnames)
  
  #select largest frac
  chr_arm_cn_fractions_highest = chr_arm_cn_fractions[
    chr_arm_cn_fractions$frac_covered == ave(chr_arm_cn_fractions$frac_covered, chr_arm_cn_fractions$patient_chr_arm,FUN=max),]
  chr_arm_cn_fractions_highest = chr_arm_cn_fractions_highest %>% dplyr::mutate(chr_arm_frac_stable = seg_covered_cr_stable > threshold_chr_arm_frac_cn_state) 
  
  chr_arm_cn_fractions_highest$region="chr_arm"
  
  #remove ones that were chrom stable
  chr_arm_cn_fractions_highest = chr_arm_cn_fractions_highest %>% 
    dplyr::mutate(selected = chr_arm_frac_stable & !patient_chrom %in% filter(chrom_cn_fractions_highest,selected)$patient_chrom) %>% as.data.frame()
  
  
  #output
  # chrom_cn_fractions_highest => call state per chrom,  select only chrom_frac_stable==T ==> selected
  # chr_arm_cn_fractions_highest => call state per chr arm, select only chr_arm_frac_stable==T, and remove those in chrom ==> selected
  # merge the dfs, rename the frac stable selection column, remove the helper id cols
  
  
  chromosome_cn_calls = rbind_no_colmatch(
    chr_arm_cn_fractions_highest %>% 
      dplyr::rename(region_frac_stable = chr_arm_frac_stable) %>% select(-patient_chr_arm,-patient_chrom),
    chrom_cn_fractions_highest %>% 
      dplyr::rename(region_frac_stable = chrom_frac_stable) %>% select(-patient_chrom))
  
  return(chromosome_cn_calls)
}


# Helper/display functions ----
get_cn_description = function(df,cna_id_col="patient_cna_id",overlap_col_id="overlap_set1_set2",amp_homloss=T,cr_amplification_threshold=2,cr_hom_loss_threshold = -1.5) {
  if(amp_homloss) { df = df %>% call_cna_amp_homloss(cr_amplification_threshold=cr_amplification_threshold,cr_hom_loss_threshold=cr_hom_loss_threshold)}
  df = df %>% summarize(cn_call = toString(call),
                        #lst_str(call) not because that sorts 
                        cn_seg_cnt = cnt_str(!!sym(cna_id_col)),
                        cn_avg = sum(cr_l2fc_50*!!sym(overlap_col_id)),
                        cn_descr = toString(paste0(call," (cr:",round(cr_l2fc_50,2),", frac:",round(!!sym(overlap_col_id),2),")")),
                        cn_call_change = grepl(",",cn_call))
  return(df)
}

call_cna_amp_homloss = function(df,cr_amplification_threshold=2,cr_hom_loss_threshold = -1.5) {
  convert_gr = F
  if(isS4(df)) { 
    df_names = names(df)
    df = as.data.frame(df) 
    convert_gr=T
  }
  df = df %>% call_cna()
  
  df = df %>% dplyr::mutate(call=ifelse(cr_l2fc_50>cr_amplification_threshold,"amp",call),
                            call=ifelse(cr_l2fc_50<(cr_hom_loss_threshold),"hom_loss",call))

  if(convert_gr) {
    df = GRanges(df)
    names(df)=df_names
  }
  return(df)
  
}

