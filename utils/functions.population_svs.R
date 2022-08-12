
## Load SV database (VCF)
# and annotate with its metadata
load_sv_database_vcf = function(sv_database_path) {
  sv_database_vcf = readVcf(sv_database_path , "hg38")
  sv_database_gr = rowRanges(sv_database_vcf)
  sv_database_gr$bp_name = names(sv_database_gr)
  sv_database_gr = sv_database_gr[order(sv_database_gr$bp_name)]
  
  
  sv_database_meta = as.data.frame(info(sv_database_vcf))
  sv_database_meta$bp_name = rownames(sv_database_meta)
  sv_database_meta = sv_database_meta[order(sv_database_meta$bp_name),]
  
  meta_tmp = as.data.frame(mcols(sv_database_gr))
  meta_tmp$bp_name = rownames(meta_tmp)
  meta_tmp = meta_tmp %>% left_join(sv_database_meta,by=c("bp_name")) 
  meta_tmp = meta_tmp[order(meta_tmp$bp_name),]
  
  mcols(sv_database_gr) = meta_tmp
  
  return(sv_database_gr)
}

make_range_sv_database = function(sv_database_gr){
  end(sv_database_gr[!is.na(sv_database_gr$END)])=sv_database_gr[!is.na(sv_database_gr$END)]$END
  return(sv_database_gr)
  
  #need to make a range from the sv_database_gr using the end coordinate
  keep_without_end = sv_database_gr[is.na(sv_database_gr$END)]
  
  sv_database_ranges = GRanges()
  for(bp_name in names(sv_database_gr[!is.na(sv_database_gr$END)])) {
    sv = sv_database_gr[bp_name]
    sv_metadata = mcols(sv)
    
    range = GRanges(seqnames = seqnames(sv),
                    IRanges(start = start(sv), end = sv_metadata$END))
    mcols(range)=sv_metadata
    names(range)=bp_name
    sv_database_ranges = c(sv_database_ranges,range)
  }
  
  sv_database_ranges = c(sv_database_ranges,keep_without_end)
  return(sv_database_ranges)
}


load_sv_population_db_overlaps = function(map_template_vars_patient, sv_databases_lst= c("nstd166","nstd186","dgv")) {
  database_overlaps = data.frame()
  for(sv_database_identifier in sv_databases_lst){
    db_map_template_vars = c(map_template_vars_patient,'${sv_database_identifier}'=sv_database_identifier)
    
    database_overlaps_path = stri_replace_all_fixed(sv_population_db_overlaps_path_template,names(db_map_template_vars), db_map_template_vars,vectorize=F)
    
    if(length(Sys.glob(database_overlaps_path))==0) { 
      print(paste0("WARNING: missing sv bp population db overlaps for patient: ",database_overlaps_path))
      next()
    }
    
    #Load overlaps 
    single_database_overlaps = read.table(database_overlaps_path,header=T,sep="\t",stringsAsFactors = F)
    single_database_overlaps = single_database_overlaps %>% dplyr::rename(patient_sv_name = set2, svtype=set2_svtype,
                                                                          sv_db_bp_name=set1, sv_db_svtype=set1_svtype)
    
    if(nrow(database_overlaps)==0){
      database_overlaps=rbind(database_overlaps,single_database_overlaps)
    } else {
      cols = names(single_database_overlaps)
      cols_collection = names(database_overlaps)
      col_diff_1 = cols_collection[!cols_collection %in% cols]
      col_diff_2 = cols[!cols %in% cols_collection]
      single_database_overlaps[,col_diff_1]=NA
      database_overlaps[,col_diff_2]=NA
      
      database_overlaps=rbind(database_overlaps,single_database_overlaps)
    }
  }
  
  return(database_overlaps)
}

filter_population_sv_db_overlaps = function(database_overlaps){
  # 50% reciprocal overlap is default, no further filtering
  # NOT strictly filter on SV type of database for 1st try 
  ## Filter out obvious mismatches like DEL/DUP 
  # minimum AF in population to be considered 0.01
  
  database_overlaps = database_overlaps %>%
    mutate(sv_db_svtype = ifelse(sv_db_svtype=="CNV" & !is.na(sv_db_variantsubtype), 
                                 ifelse(sv_db_variantsubtype %in% c("deletion","loss"), "DEL", 
                                        ifelse(sv_db_variantsubtype %in% c("duplication","gain"), "DUP",sv_db_svtype)),sv_db_svtype))
  
  #helps to reduce ambiguous cnv count
  #database_overlaps %>% filter(sv_db_svtype=="CNV") %>% nrow()
  
  ## Filter out obvious mismatches like DEL/DUP 
  database_overlaps = database_overlaps %>% filter( !( (sv_db_svtype=="DEL" & svtype=="DUP") |  (sv_db_svtype=="DUP" & svtype=="DEL")))
  
  
  #make AF for dgv
  database_overlaps = database_overlaps %>% 
    mutate(sv_db_af = ifelse(is.na(sv_db_af), (observedgains+observedlosses)/samplesize, sv_db_af)) %>%
    filter(!is.na(sv_db_af) & !(sv_db_af==1 & samplesize<50))
  #database_overlaps %>% filter(observedgains>10 & observedlosses>10) %>% select(sv_db_variantsubtype) %>% unique()
  #sv_db_variantsubtype is gain+loss or  complex
  
  if(FALSE) {
    
    #potential rescue in case overlaps multiple low AF variants
    #summed_population_af = database_overlaps %>% filter(sv_db_af<0.01) %>% 
    # select(patient_sv_name,sv_db_bp_name,sv_database,sv_db_af) %>% unique() %>% 
    # group_by(patient_sv_name,sv_database) %>% summarize(sv_db_af_summed=sum(sv_db_af)) %>%
    # filter(sv_db_af_summed>0.01)
  }
  
  # minimum AF in population to be considered 0.01
  #database_overlaps_bk=database_overlaps
  database_overlaps = database_overlaps %>% filter(sv_db_af>=0.01)  
  
  if(FALSE) {
    #remove the ones that already overlapped a good one
    #summed_population_af = summed_population_af %>% filter(!patient_sv_name %in% database_overlaps$patient_sv_name)
    #what would be rescued and do I agree?
    #database_overlaps_bk %>% filter(patient_sv_name %in% summed_population_af$patient_sv_name) %>% arrange(patient_sv_name) %>% head()
    #do not think it is a good idea. Seems that variants can be in databases  multiple times eg.g as insertion or duplication
    
    #svs backup for comparison
    #svs_df_no_overlaps = svs_df %>% filter(!patient_sv_name %in% database_overlaps$patient_sv_name)
    
  }
  
  
  return(database_overlaps)
}
