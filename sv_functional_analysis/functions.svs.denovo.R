

read_delly_sv_vcf_trio = function(vcf_path,trio) {
  allele_fraction_cols = c("proband_af","father_af","mother_af")
  metadata_cols = c("sourceId",allele_fraction_cols)
  
  vcf = readVcf(vcf_path, "hg38") 
  
  if(length(vcf)==0){ return(GRanges()) }
  
  #prevent error
  vcf = vcf[!isNA(info(vcf)$SVTYPE)]
  info(vcf)[info(vcf)$SVTYPE=="BND",c("SVTYPE")]="TRA"
  
  ## Make dataframe based on this
  all_gr = breakpointRanges(vcf) 
  
  ## sourceId for row names, which splits bp_nm into _bp1 and _bp2
  gr_df = as.data.frame(all_gr) %>% dplyr::mutate(bp_name = rownames(.)) %>%  dplyr::select(sourceId)
  
  # adding names from row ranges because also not filtered yet, and often missing if subsetted vcf
  vcf_geno = as.data.frame(info(vcf)) 
  vcf_geno$sourceId = names(rowRanges(vcf))
  vcf_geno  = vcf_geno %>%  dplyr::select(sourceId, PRECISE)
  
  #filter based on if occurs in germline_gr => some unpaired were removed during loading
  vcf_geno = vcf_geno %>% filter(sourceId %in% gr_df$sourceId)
  
  ## Calculate AF: how to do this depends on precise/imprecise, which is in info tag 
  #AF = RV/(RR+RV) for precise variants
  #AF = DV/(DR+DV) for imprecise variants
  vcf_geno_precise = 
    filter(vcf_geno,PRECISE) %>% 
    left_join (
      as.data.frame(geno(vcf)$RV) %>% dplyr::mutate(sourceId = rownames(.)) %>% 
        dplyr::rename("proband_RV" := !!trio$proband_id) %>% dplyr::rename("father_RV" := !!trio$father_id) %>% dplyr::rename("mother_RV" := !!trio$mother_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(vcf)$RR) %>% dplyr::mutate(sourceId = rownames(.)) %>%
        dplyr::rename("proband_RR" := !!trio$proband_id) %>% dplyr::rename("father_RR" := !!trio$father_id) %>% dplyr::rename("mother_RR" := !!trio$mother_id)
      , by="sourceId"
    ) %>%
    dplyr::mutate(proband_af = (proband_RV/(proband_RR+proband_RV))) %>%
    dplyr::mutate(father_af = (father_RV/(father_RR+father_RV))) %>%
    dplyr::mutate(mother_af = (mother_RV/(mother_RR+mother_RV)))
  
  #if both are 0 then drop the NAs, they can be rescued by DVs or otherwise we just dont have data for that variant
  vcf_geno_precise = drop_na(vcf_geno_precise)
  
  # calculate based on DV for all to also rescue the NAs with RV
  vcf_geno_all = 
    vcf_geno %>% 
    left_join (
      as.data.frame(geno(vcf)$DV) %>% dplyr::mutate(sourceId = rownames(.)) %>% 
        dplyr::rename("proband_DV" := !!trio$proband_id) %>% dplyr::rename("father_DV" := !!trio$father_id) %>% dplyr::rename("mother_DV" := !!trio$mother_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(vcf)$DR) %>% dplyr::mutate(sourceId = rownames(.)) %>%
        dplyr::rename("proband_DR" := !!trio$proband_id) %>% dplyr::rename("father_DR" := !!trio$father_id) %>% dplyr::rename("mother_DR" := !!trio$mother_id)
      , by="sourceId"
    ) %>%
    dplyr::mutate(proband_af = (proband_DV/(proband_DR+proband_DV))) %>%
    dplyr::mutate(father_af = (father_DV/(father_DR+father_DV))) %>%
    dplyr::mutate(mother_af = (mother_DV/(mother_DR+mother_DV)))
  
  vcf_geno_all = drop_na(vcf_geno_all)
  
  #if AF in precise table, then exclude it here 
  vcf_af = rbind(vcf_geno_precise[,metadata_cols], 
                 vcf_geno_all[!vcf_geno_all$sourceId %in% vcf_geno_precise$sourceId,metadata_cols])
  
  read_support_colnames =  names(vcf_geno_all)[grepl("DV|DR|RV|RR",names(vcf_geno_all))]
  read_support_colnames_precise =  names(vcf_geno_precise)[grepl("DV|DR|RV|RR",names(vcf_geno_precise))]
  
  vcf_af = vcf_af %>% left_join(vcf_geno_all[,c("sourceId",read_support_colnames)],by="sourceId") %>%
    left_join(vcf_geno_precise[,c("sourceId",read_support_colnames_precise)],by="sourceId") 
  
  #merge allele frequencies, 
  gr_df = gr_df %>% left_join(vcf_af, by="sourceId")
  
  #add the annotation to breakpoints
  read_support_colnames =  names(vcf_af)[grepl("DV|DR|RV|RR",names(vcf_af))]
  
  all_gr = annotate_metadata(all_gr,gr_df,metadata_cols = unique(c(read_support_colnames,allele_fraction_cols)))
  
  #annotate with tool
  mcols(all_gr)[["tool"]]="delly"
  
  return(all_gr)
}
read_gridss_sv_vcf_trio = function(vcf_path,trio) {
  allele_fraction_cols = c("proband_af","father_af","mother_af")
  
  vcf = readVcf(vcf_path, "hg38")
  if(length(vcf)==0){ return(GRanges()) }
  
  all_gr = breakpointRanges(vcf)
  elementMetadata(all_gr)["svtype"] = get_svtype(all_gr)
  
  #this removes unpartnered, need to subset afterwards
  vcf = vcf[names(all_gr)]
  
  vcf_geno = geno(vcf)
  gr_df = as.data.frame(all_gr) %>% dplyr::mutate(bp_name = rownames(.))
  
  ## Allele frequency calculation
  #how to calculate AF depends on variant size 
  # svlen <1000 dont use the REFPAIR otherwise do 
  #uses the label of patient not the id.
  
  vcf_anno = gr_df %>% dplyr::select(bp_name,svLen) %>% left_join(
    as.data.frame(vcf_geno$VF) %>% dplyr::mutate(bp_name = rownames(.)) %>% 
      dplyr::rename("proband_VF" := !!trio$proband_label) %>% dplyr::rename("father_VF" := !!trio$father_label) %>% dplyr::rename("mother_VF" := !!trio$mother_label)
    , by="bp_name") %>% left_join(
      as.data.frame(vcf_geno$REF) %>% dplyr::mutate(bp_name = rownames(.)) %>% 
        dplyr::rename("proband_REF" := !!trio$proband_label) %>% dplyr::rename("father_REF" := !!trio$father_label) %>% dplyr::rename("mother_REF" := !!trio$mother_label)
      , by="bp_name") %>% left_join(
        as.data.frame(vcf_geno$REFPAIR) %>% dplyr::mutate(bp_name = rownames(.)) %>% 
          dplyr::rename("proband_REFPAIR" := !!trio$proband_label) %>% dplyr::rename("father_REFPAIR" := !!trio$father_label) %>% dplyr::rename("mother_REFPAIR" := !!trio$mother_label)
        , by="bp_name") %>%
    dplyr::mutate(proband_af = ifelse(!is.na(svLen)&svLen<1000, (proband_VF/(proband_VF+proband_REF)), (proband_VF/(proband_VF+proband_REF+proband_REFPAIR)))) %>%
    dplyr::mutate(father_af = ifelse(!is.na(svLen)&svLen<1000, (father_VF/(father_VF+father_REF)), (father_VF/(father_VF+father_REF+father_REFPAIR))))  %>%
    dplyr::mutate(mother_af = ifelse(!is.na(svLen)&svLen<1000, (mother_VF/(mother_VF+mother_REF)), (mother_VF/(mother_VF+mother_REF+mother_REFPAIR))))
  
  vcf_anno[is.na(vcf_anno)]=0
  #add the annotation to breakpoints
  read_support_colnames =  names(vcf_anno)[grepl("VF|REF",names(vcf_anno))]
  
  all_gr = annotate_metadata(all_gr,vcf_anno,metadata_cols = c(read_support_colnames,allele_fraction_cols))
  
  #annotate with tool
  mcols(all_gr)[["tool"]]="gridss"
  return(all_gr)
}

manta_metadata_trio = function(vcf, support_type, trio) {
  
  sample_ids=c(trio$proband_id,trio$father_id,trio$mother_id)
  
  vcf_geno_df = as.data.frame(geno(vcf)[[support_type]])
  vcf_geno_df$sourceId = rownames(vcf_geno_df)
  
  vcf_geno_df[,sample_ids] = 
    lapply(vcf_geno_df[,sample_ids], function(x){gsub("[c()]", "", x)})
  
  vcf_geno_df[,sample_ids] = 
    lapply(vcf_geno_df[,sample_ids], function(x){gsub(":", ", ",  x)})
  
  vcf_geno_df[vcf_geno_df == "integer0"]=NA
  
  
  column_label="proband"
  vcf_geno_df = separate(vcf_geno_df, col =trio$proband_id, 
                         into = c(paste0(column_label,"_REF_",support_type),paste0(column_label,"_VAR_",support_type)), sep = ", ")
  
  column_label="father"
  vcf_geno_df = separate(vcf_geno_df, col =trio$father_id, 
                         into = c(paste0(column_label,"_REF_",support_type),paste0(column_label,"_VAR_",support_type)), sep = ", ")
  
  column_label="mother"
  vcf_geno_df = separate(vcf_geno_df, col =trio$mother_id, 
                         into = c(paste0(column_label,"_REF_",support_type),paste0(column_label,"_VAR_",support_type)), sep = ", ")
  
  
  
  return(vcf_geno_df)
}  

read_manta_sv_vcf_trio = function(vcf_path, trio) {
  allele_fraction_cols = c("proband_af","father_af","mother_af")
  
  vcf = readVcf(vcf_path, "hg38")
  vcf_geno=data.frame()
  
  if(length(vcf)==0){
    return(Granges)
  }
  
  #SR and PR separate dataframes at first
  vcf_geno_PR = manta_metadata_trio(vcf, "PR", trio)
  vcf_geno_SR = manta_metadata_trio(vcf, "SR", trio)
  
  ## join and then use both for calculation of AF 
  vcf_geno = vcf_geno_PR %>% left_join(vcf_geno_SR,by="sourceId") 
  
  read_support_colnames =  names(vcf_geno)[grepl("PR|SR",names(vcf_geno))]
  
  vcf_geno[is.na(vcf_geno)]=0
  vcf_geno[,read_support_colnames] = as.data.frame(sapply(vcf_geno[,read_support_colnames], as.numeric))
  
  vcf_geno$proband_af = (vcf_geno$proband_VAR_PR + vcf_geno$proband_VAR_SR) / 
    (vcf_geno$proband_REF_SR + vcf_geno$proband_REF_PR +vcf_geno$proband_VAR_SR + vcf_geno$proband_VAR_PR) 
  
  vcf_geno$father_af = (vcf_geno$father_VAR_PR + vcf_geno$father_VAR_SR) / 
    (vcf_geno$father_REF_SR + vcf_geno$father_REF_PR +vcf_geno$father_VAR_SR + vcf_geno$father_VAR_PR) 
  
  vcf_geno$mother_af = (vcf_geno$mother_VAR_PR + vcf_geno$mother_VAR_SR) / 
    (vcf_geno$mother_REF_SR + vcf_geno$mother_REF_PR +vcf_geno$mother_VAR_SR + vcf_geno$mother_VAR_PR) 
  
  # in case the source Id and bp name do not match then the somatic VCF geno "bp name" matches the sourceId
  ## vcf_geno[!vcf_geno$sourceId %in% gr$sourceId & !vcf_geno$sourceId %in% names(gr),]
  ### sometimes another does match with the same info like MantaDUP:TANDEM:35919:0:1:0:0:0 :1/2/3 etc
  ##vcf_geno[vcf_geno$sourceId %in% (gr[names(gr) != gr$sourceId ])$sourceId,]
  
  vcf_geno[is.na(vcf_geno)]=0
  
  
  vcf_geno$sourceId= gsub(".",":",vcf_geno$sourceId,fixed=T)
  
  
  all_gr = breakpointRanges(vcf)
  all_gr = all_gr[grepl("Manta",all_gr$sourceId),]
  all_gr = annotate_metadata(all_gr,vcf_geno,metadata_cols = c(read_support_colnames,allele_fraction_cols))
  
  #add another column for analysis type
  #annotate with tool
  if(length(all_gr)>0)  mcols(all_gr)[["tool"]]="manta"
  return(all_gr) 
}



annotate_sv_af_class_trio = function(df,threshold_presence=0.1, threshold_absense=0.05) {
  svs_df$variant_type = NULL
  #proband only #denovo
  #father only
  #mother only 
  #proband father
  #proband mother
  #shared 
  
  df = df %>% dplyr::mutate(variant_type=
                              ifelse(  proband_af>=threshold_presence & father_af>=threshold_presence & mother_af>=threshold_presence ,"shared",
                                       ifelse( proband_af>=threshold_presence & father_af>=threshold_presence  & mother_af<threshold_absense,"proband_father",
                                               ifelse( proband_af>=threshold_presence & mother_af>=threshold_presence & father_af<threshold_absense,"proband_mother",
                                                       ifelse( mother_af>=threshold_presence & proband_af<threshold_absense & father_af<threshold_absense,"mother_only",
                                                               ifelse( father_af>=threshold_presence & proband_af<threshold_absense & mother_af<threshold_absense,"father_only",
                                                                       ifelse( proband_af>=threshold_presence & father_af<threshold_absense & mother_af<threshold_absense,"proband_only",
                                                                               ifelse( proband_af<threshold_presence & father_af<threshold_absense & mother_af<threshold_absense,"proband_low_af",
                                                                                       "multiple_low_af"))))))))
  df$variant_type = factor(df$variant_type)
  
  return(df)
}
