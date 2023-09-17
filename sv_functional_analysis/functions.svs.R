## VCF loading functions
get_svtype <- function(gr) {
  #Function copied from GRIDSS example code
  # CTX if seq names are the same
  ## return as complex if unpartnered
  
  return_gr = gr
  unpartnered_gr=gr[!gr$partner %in% names(gr)]
  gr = gr[gr$partner %in% names(gr)]
  
  svtype = ifelse(seqnames(gr) != seqnames(gr[gr$partner]), "CTX", # inter-chromosomosal
                  ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                         ifelse(strand(gr) == strand(gr[gr$partner]), "INV",
                                ifelse(xor(start(gr) < start(gr[gr$partner]), strand(gr) == "-"), "DEL",
                                       "DUP"))))
  
  if(length(unpartnered_gr)>0){
    return_gr[names(unpartnered_gr)]$svtype="complex"
  }
  return_gr[names(gr)]$svtype=svtype
  return(  return_gr$svtype)
}


get_ctx_strands = function(merged_svs) {
  merged_svs = merged_svs %>% rowwise() %>% 
    mutate(strand_source_fw= (svtype=="CTX"&grepl("+",sv_merged_coordinate,fixed = T)), 
           strand_partner_fw= (svtype=="CTX"&grepl("+",partner_sv_merged_coordinate,fixed = T)), 
           ctx_strands = ifelse(svtype!="CTX",NA,
                                ifelse(strand_source_fw&strand_partner_fw,"+/+",
                                       ifelse(strand_source_fw&!strand_partner_fw,"+/-",
                                              ifelse(!strand_source_fw&strand_partner_fw,"-/+","-/-")))))
  
  return(merged_svs)
}


get_ctx_types = function(merged_svs,group_cols=c("patient_label","complex_sv_id")) {
  if(!"ctx_strands" %in% names(merged_svs)) {
    print("need ctx_strands")
    return()
  }
  ctx_types = merged_svs %>% 
    group_by(across(all_of(group_cols))) %>%
    summarize(ctx_inv_cnt=sum(ctx_strands=="+/+" | ctx_strands=="-/-"),
              ctx_dup_cnt=sum(ctx_strands=="-/+"),
              ctx_del_cnt=sum(ctx_strands=="+/-"))
  
  return(ctx_types)
}


## copied from match wgs functions but made safe for ensdb
annotate_metadata = function(all_gr,metadata_df,metadata_cols=c("tumor_af","normal_af","somatic")) {
  all_gr=all_gr[order(names(all_gr))]
  metadata = mcols(all_gr,use.names = T)  
  metadata$bp_name = rownames(metadata)
  
  if("sourceId" %in% names(metadata_df) && nrow(dplyr::filter(metadata_df,sourceId %in% metadata$sourceId))>0) {
    metadata = metadata %>% merge(metadata_df %>% dplyr::select(all_of(c("sourceId",metadata_cols))),by="sourceId",all.x=T) %>% unique()
  } 
  else if("bp_name" %in% names(metadata_df) && nrow(dplyr::filter(metadata_df,bp_name %in% metadata$bp_name))>0) {
    metadata = metadata %>% merge(metadata_df %>% dplyr::select(all_of(c("bp_name",metadata_cols))),by="bp_name",all.x=T) %>% unique()
  }
  metadata=metadata[order(metadata$bp_name),]
  rownames(metadata) = metadata$bp_name  
  
  mcols(all_gr) <- metadata
  return(all_gr)
}
manta_metadata = function(vcf, patient, support_type, somatic=TRUE) {
  #vcf object as input
  vcf_geno_df = as.data.frame(geno(vcf)[[support_type]])
  vcf_geno_df$sourceId = rownames(vcf_geno_df)
  
  vcf_geno_df[,as.character(patient$normal_id)] = 
    gsub("[c()]", "", vcf_geno_df[,as.character(patient$normal_id)])
  
  vcf_geno_df[,as.character(patient$normal_id)] = 
    gsub(":", ", ", vcf_geno_df[,as.character(patient$normal_id)])
  
  vcf_geno_df[vcf_geno_df == "integer0"]=NA
  
  vcf_geno_df = separate(vcf_geno_df, col = as.character(patient$normal_id), 
                         into = c(paste0("normal_REF_",support_type),paste0("normal_VAR_",support_type)), sep = ", ")
  
  #remove the ":0" if couldnt be split 
  #vcf_geno_df[,paste0("normal_REF_",support_type)] = sub(":0", "", vcf_geno_df[,paste0("normal_REF_",support_type)])
  
  if(somatic){
    vcf_geno_df[,as.character(patient$tumor_id)] = 
      gsub("[c()]", "", vcf_geno_df[,as.character(patient$tumor_id)])
    
    vcf_geno_df[,as.character(patient$tumor_id)] = 
      gsub(":", ", ", vcf_geno_df[,as.character(patient$tumor_id)])
    
    vcf_geno_df[vcf_geno_df == "integer0"]=NA
    
    vcf_geno_df = separate(vcf_geno_df, col = as.character(patient$tumor_id), 
                           into = c(paste0("tumor_REF_",support_type),paste0("tumor_VAR_",support_type)), sep = ", ")
    
    #remove the ":0" if couldnt be split 
    # vcf_geno_df[,paste0("tumor_REF_",support_type)] = sub(":0", "", vcf_geno_df[,paste0("tumor_REF_",support_type)])
  }
  
  
  #vcf_geno_df[is.na(vcf_geno_df)]=0
  
  return(vcf_geno_df)
}  

read_manta_sv_vcf = function(vcf_germline_path, vcf_somatic_path,patient) {
  
  # SOMATIC and GERMLINE are separate files for somatic and germline breakpoints
  ## join and then use both for calculation of AF 
  
  somatic_vcf = readVcf(vcf_somatic_path, "hg38")
  somatic_vcf_geno=data.frame()
  
  if(length(somatic_vcf)>0){ 
    
    #SR and PR separate dataframes at first
    somatic_vcf_geno_PR = manta_metadata(somatic_vcf, patient, "PR", somatic=TRUE)
    somatic_vcf_geno_SR = manta_metadata(somatic_vcf, patient, "SR", somatic=TRUE)
    
    ## join and then use both for calculation of AF 
    somatic_vcf_geno = somatic_vcf_geno_SR %>% left_join(somatic_vcf_geno_PR,by="sourceId")
    
    somatic_vcf_geno$somatic=TRUE
    
    read_support_colnames =  names(somatic_vcf_geno)[grepl("PR|SR",names(somatic_vcf_geno))]
    somatic_vcf_geno[is.na(somatic_vcf_geno)]=0
    somatic_vcf_geno[,read_support_colnames] = as.data.frame(sapply(somatic_vcf_geno[,read_support_colnames], as.numeric))

    somatic_vcf_geno$tumor_af = (somatic_vcf_geno$tumor_VAR_PR + somatic_vcf_geno$tumor_VAR_SR) / 
      (somatic_vcf_geno$tumor_REF_SR + somatic_vcf_geno$tumor_REF_PR +somatic_vcf_geno$tumor_VAR_SR + somatic_vcf_geno$tumor_VAR_PR) 
    
    somatic_vcf_geno$normal_af = (somatic_vcf_geno$normal_VAR_PR + somatic_vcf_geno$normal_VAR_SR) / 
      (somatic_vcf_geno$normal_REF_SR + somatic_vcf_geno$normal_REF_PR +somatic_vcf_geno$normal_VAR_SR + somatic_vcf_geno$normal_VAR_PR) 
    
    # in case the source Id and bp name do not match then the somatic VCF geno "bp name" matches the sourceId
    ## somatic_vcf_geno[!somatic_vcf_geno$sourceId %in% somatic_gr$sourceId & !somatic_vcf_geno$sourceId %in% names(somatic_gr),]
    ### sometimes another does match with the same info like MantaDUP:TANDEM:35919:0:1:0:0:0 :1/2/3 etc
    ##somatic_vcf_geno[somatic_vcf_geno$sourceId %in% (somatic_gr[names(somatic_gr) != somatic_gr$sourceId ])$sourceId,]
    
    somatic_vcf_geno$sourceId= gsub(".",":",somatic_vcf_geno$sourceId,fixed=T)
    
  } 
  
  ## germline has only normal sample 
  germline_vcf = readVcf(vcf_germline_path, "hg38")
  germline_vcf_geno=data.frame()
  
  if(length(germline_vcf)>0) {
    
    germline_vcf_geno_PR = manta_metadata(germline_vcf, patient, "PR", somatic=FALSE)
    germline_vcf_geno_SR = manta_metadata(germline_vcf, patient, "SR", somatic=FALSE)
    
    germline_vcf_geno = germline_vcf_geno_SR %>% left_join(germline_vcf_geno_PR,by="sourceId")
    germline_vcf_geno$somatic=FALSE
    
    ## lacking tumor AF so adding dummy variable for merging
    germline_vcf_geno$tumor_af = NA
    
    read_support_colnames =  names(germline_vcf_geno)[grepl("PR|SR",names(germline_vcf_geno))]
    germline_vcf_geno[is.na(germline_vcf_geno)]=0
    germline_vcf_geno[,read_support_colnames] = as.data.frame(sapply(germline_vcf_geno[,read_support_colnames], as.numeric))
    
    germline_vcf_geno$normal_af = (germline_vcf_geno$normal_VAR_PR + germline_vcf_geno$normal_VAR_SR) / 
      (germline_vcf_geno$normal_REF_SR + germline_vcf_geno$normal_REF_PR +germline_vcf_geno$normal_VAR_SR + germline_vcf_geno$normal_VAR_PR) 
    
    germline_vcf_geno$sourceId= gsub(".",":",germline_vcf_geno$sourceId,fixed=T)
    
  }
  
  ## Note: manta sometimes 0s for normal AF in diploid model, also somtimes normal > tumor AF in somatic 
  
  #breakpoints need to append and annotate with allele frequencies, tool, somatic or not 
  ##  filtering the sourceIds for "manta" to remove "svrecord" since they ofen overlap with normal variants, but without metadata
  
  
  read_support_colnames_somatic =  names(somatic_vcf_geno)[grepl("PR|SR",names(somatic_vcf_geno))]
  
  somatic_gr = breakpointRanges(somatic_vcf)
  somatic_gr = somatic_gr[grepl("Manta",somatic_gr$sourceId),]
  somatic_gr = annotate_metadata(somatic_gr,somatic_vcf_geno,metadata_cols = c(read_support_colnames_somatic,"tumor_af","normal_af","somatic"))
  
  germline_gr = breakpointRanges(germline_vcf)
  germline_gr = germline_gr[grepl("Manta",germline_gr$sourceId),]
  
  read_support_colnames_germline =  names(germline_vcf_geno)[grepl("PR|SR",names(germline_vcf_geno))]
  
  missing_cols=read_support_colnames_somatic[!read_support_colnames_somatic %in% read_support_colnames_germline]
  germline_vcf_geno[,missing_cols]=NA
  
  germline_gr = annotate_metadata(germline_gr,germline_vcf_geno,metadata_cols = c(read_support_colnames_somatic,"tumor_af","normal_af","somatic"))
    
  all_gr = c(somatic_gr,germline_gr)
  #add another column for analysis type
  #annotate with tool
  if(length(all_gr)>0)  mcols(all_gr)[["tool"]]="manta"
  return(all_gr) 
}
read_delly_sv_vcf = function(vcf_path,patient) {
  ##NOTE: somatic disabled because  DV/DR/RV/RR and allele frequencies are these same
  ## Function is available for tagging based on name but not currently used by Fusion-sq 
  
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
        dplyr::rename("tumor_RV" := !!patient$tumor_id) %>% dplyr::rename("normal_RV" := !!patient$normal_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(vcf)$RR) %>% dplyr::mutate(sourceId = rownames(.)) %>%
        dplyr::rename("tumor_RR" := !!patient$tumor_id) %>% dplyr::rename("normal_RR" := !!patient$normal_id)
      , by="sourceId"
    ) %>%
    dplyr::mutate(tumor_af = (tumor_RV/(tumor_RR+tumor_RV))) %>%
    dplyr::mutate(normal_af = (normal_RV/(normal_RR+normal_RV)))
  
  #if both are 0 then drop the NAs, they can be rescued by DVs or otherwise we just dont have data for that variant
  vcf_geno_precise = drop_na(vcf_geno_precise)
  
  # calculate based on DV for all to also rescue the NAs with RV
  vcf_geno_all = 
    vcf_geno %>% 
    left_join (
      as.data.frame(geno(vcf)$DV) %>% dplyr::mutate(sourceId = rownames(.)) %>% 
        dplyr::rename("tumor_DV" := !!patient$tumor_id) %>% dplyr::rename("normal_DV" := !!patient$normal_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(vcf)$DR) %>% dplyr::mutate(sourceId = rownames(.)) %>%
        dplyr::rename("tumor_DR" := !!patient$tumor_id) %>% dplyr::rename("normal_DR" := !!patient$normal_id)
      , by="sourceId"
    ) %>%
    dplyr::mutate(tumor_af = (tumor_DV/(tumor_DR+tumor_DV))) %>%
    dplyr::mutate(normal_af = (normal_DV/(normal_DR+normal_DV)))
  
  vcf_geno_all = drop_na(vcf_geno_all)
  
  metadata_cols = c("sourceId","tumor_af","normal_af")
  #if AF in precise table, then exclude it here 
  vcf_af = rbind(vcf_geno_precise[,metadata_cols], 
                 vcf_geno_all[!vcf_geno_all$sourceId %in% vcf_geno_precise$sourceId,metadata_cols])
  
  read_support_colnames =  names(vcf_geno_all)[grepl("DV|DR|RV|RR",names(vcf_geno_all))]
  read_support_colnames_precise =  names(vcf_geno_precise)[grepl("DV|DR|RV|RR",names(vcf_geno_precise))]
  
  vcf_af = vcf_af %>% left_join(vcf_geno_all[,c("sourceId",read_support_colnames)],by="sourceId") %>%
    left_join(vcf_geno_precise[,c("sourceId",read_support_colnames_precise)],by="sourceId") 
  
  ## DISABLED: If somatic file is provided, flag breakpoints as such by annotate metadata
  #input germline anno table, somatic file path, total intervals = GRanges() for delly 
  ##use the full DF for annotation here instead of the allele frequency one because of "drop NA"
  #gr_df = annotate_metadata_somatic(gr_df, argv$vcf_somatic, GRanges())
  
  ## Currently not used:
  gr_df$somatic = NA
  
  #merge allele frequencies, 
  gr_df = gr_df %>% left_join(vcf_af, by="sourceId")
  
  #add the annotation to breakpoints
  read_support_colnames =  names(vcf_af)[grepl("DV|DR|RV|RR",names(vcf_af))]
  
  all_gr = annotate_metadata(all_gr,gr_df,metadata_cols = c(read_support_colnames,"tumor_af","normal_af","somatic"))
  
  
  #annotate with tool
  mcols(all_gr)[["tool"]]="delly"
  
  return(all_gr)
}
read_gridss_sv_vcf = function(vcf_path,patient) {
  ##NOTE: somatic disabled because allele frequencies are these same
  ## Function is available for tagging based on name but not currently used by Fusion-sq 
  
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
  #uses the label of tumor/normal not the id.
  
  vcf_anno = gr_df %>% dplyr::select(bp_name,svLen) %>% left_join(
    as.data.frame(vcf_geno$VF) %>% dplyr::mutate(bp_name = rownames(.)) %>% 
      dplyr::rename("tumor_VF" := !!patient$tumor_label) %>% dplyr::rename("normal_VF" := !!patient$normal_label)
    , by="bp_name") %>% left_join(
      as.data.frame(vcf_geno$REF) %>% dplyr::mutate(bp_name = rownames(.)) %>% 
        dplyr::rename("tumor_REF" := !!patient$tumor_label) %>% dplyr::rename("normal_REF" := !!patient$normal_label)
      , by="bp_name") %>% left_join(
        as.data.frame(vcf_geno$REFPAIR) %>% dplyr::mutate(bp_name = rownames(.)) %>% 
          dplyr::rename("tumor_REFPAIR" := !!patient$tumor_label) %>% dplyr::rename("normal_REFPAIR" := !!patient$normal_label)
        , by="bp_name") %>%
    dplyr::mutate(tumor_af = ifelse(!is.na(svLen)&svLen<1000, (tumor_VF/(tumor_VF+tumor_REF)), (tumor_VF/(tumor_VF+tumor_REF+tumor_REFPAIR)))) %>%
    dplyr::mutate(normal_af = ifelse(!is.na(svLen)&svLen<1000, (normal_VF/(normal_VF+normal_REF)), (normal_VF/(normal_VF+normal_REF+normal_REFPAIR))))
  
  ## DISABLED: annotate with somatic flag
  #vcf_anno = annotate_metadata_somatic(vcf_anno, argv$vcf_somatic, total_intervals)
  #rownames(vcf_anno) = vcf_anno$bp_name
  vcf_anno$somatic=NA
  
  
  #add the annotation to breakpoints
  read_support_colnames =  names(vcf_anno)[grepl("VF|REF",names(vcf_anno))]
  
  all_gr = annotate_metadata(all_gr,vcf_anno,metadata_cols = c(read_support_colnames,"tumor_af","normal_af","somatic"))

  #annotate with tool
  mcols(all_gr)[["tool"]]="gridss"
  return(all_gr)
}

#germline only analysis, nienkes project
## DEPRECATED use read_manta_diploidSV()
load_manta_diploidSV = function(vcf_path,sample) {
  
  germline_vcf = readVcf(vcf_path, "hg38")
  
  
  ## join and then use both for calculation of AF 
  germline_vcf_geno_PR = manta_metadata(germline_vcf, sample, "PR", somatic=FALSE)
  germline_vcf_geno_SR = manta_metadata(germline_vcf, sample, "SR", somatic=FALSE)
  
  germline_vcf_geno = germline_vcf_geno_SR %>% left_join(germline_vcf_geno_PR,by="sourceId")
  
  ## lacking tumor AF so adding dummy variable for merging
  germline_vcf_geno$tumor_af = NA
  germline_vcf_geno$somatic = NA
  germline_vcf_geno$normal_af = (as.integer(germline_vcf_geno$normal_VAR_PR) + as.integer(germline_vcf_geno$normal_VAR_SR)) / 
    ( as.integer(germline_vcf_geno$normal_REF_SR) + as.integer(germline_vcf_geno$normal_REF_PR) +
        as.integer(germline_vcf_geno$normal_VAR_SR) + as.integer(germline_vcf_geno$normal_VAR_PR) )
  
  #breakpoints need to append and annotate with allele frequencies, tool, somatic or not 
  germline_gr = breakpointRanges(germline_vcf)
  
  ## Note: manta sometimes 0s for normal AF in diploid model, also somtimes normal > tumor AF in somatic 
  
  ##  filtering the sourceIds for "manta" to remove "svrecord" since they ofen overlap with normal variants, but without metadata
  
  # in case the source Id and bp name do not match then the somatic VCF geno "bp name" matches the sourceId
  germline_gr = germline_gr[grepl("Manta",germline_gr$sourceId),]
  germline_vcf_geno$sourceId= gsub(".",":",germline_vcf_geno$sourceId,fixed=T)
  germline_gr = annotate_metadata(germline_gr,germline_vcf_geno)
  return(germline_gr)
}


read_manta_diploidSV = function(vcf_germline_path,patient) {
  ## based on read_manta_sv_vcf() so fixes the tumor/normal AF 
  
  
  germline_vcf = readVcf(vcf_germline_path, "hg38")
  germline_vcf_geno=data.frame()
  
  if(length(germline_vcf)==0) {
    return(GRanges())
  }
    germline_vcf_geno_PR = manta_metadata(germline_vcf, patient, "PR", somatic=FALSE)
    germline_vcf_geno_SR = manta_metadata(germline_vcf, patient, "SR", somatic=FALSE)
    
    germline_vcf_geno = germline_vcf_geno_SR %>% left_join(germline_vcf_geno_PR,by="sourceId")
    germline_vcf_geno$somatic=FALSE
    
    ## lacking tumor AF so adding dummy variable for merging
    germline_vcf_geno$tumor_af = NA
    
    read_support_colnames =  names(germline_vcf_geno)[grepl("PR|SR",names(germline_vcf_geno))]
    germline_vcf_geno[is.na(germline_vcf_geno)]=0
    germline_vcf_geno[,read_support_colnames] = as.data.frame(sapply(germline_vcf_geno[,read_support_colnames], as.numeric))
    
    germline_vcf_geno$normal_af = (germline_vcf_geno$normal_VAR_PR + germline_vcf_geno$normal_VAR_SR) / 
      (germline_vcf_geno$normal_REF_SR + germline_vcf_geno$normal_REF_PR +germline_vcf_geno$normal_VAR_SR + germline_vcf_geno$normal_VAR_PR) 
    
    germline_vcf_geno$sourceId= gsub(".",":",germline_vcf_geno$sourceId,fixed=T)
    
  
  
 
  germline_gr = breakpointRanges(germline_vcf)
  germline_gr = germline_gr[grepl("Manta",germline_gr$sourceId),]
  
  germline_gr = annotate_metadata(germline_gr,germline_vcf_geno,metadata_cols = c(read_support_colnames,"tumor_af","normal_af","somatic"))
  
  
  if(length(germline_gr)>0)  mcols(germline_gr)[["tool"]]="manta"
  return(germline_gr) 
}
### ENDOF VCF loading functions

## Variant classification
#deprecated? 2021-08-19  => see sep functions file
annotate_variant_af_class = function(summary_df) {
  summary_df = summary_df %>% 
    dplyr::mutate(tumor_normal_diff = (tumor_af-normal_af)) %>%
    dplyr::mutate(tumor_normal_ratio = (tumor_normal_diff/normal_af)) 
  
  summary_df = summary_df %>% dplyr::mutate(somatic_variant = (tumor_normal_diff>0.05 & tumor_normal_ratio>1.5))
  summary_df = summary_df %>% dplyr::mutate(germline_variant = (normal_af>0.05 & tumor_normal_ratio<1.1))
  summary_df = summary_df %>% dplyr::mutate(low_af = (tumor_normal_diff<0.05 & normal_af<0.05))
  
  summary_df = summary_df %>% dplyr::mutate(somatic_variant = ifelse( (is.na(tumor_af)&!is.na(normal_af)), FALSE,somatic_variant))
  
  summary_df = summary_df %>% dplyr::select(-tumor_normal_diff)
  return(summary_df)
}

## End of variant classification


## Functions considering SVs as ranges, to determine same-ness of SVs and overlaps with features
df_to_gr = function(set) {
  if(is.data.frame(set)) { 
    set = unique(set)
    row.names(set) = set$bp_name
    set =GRanges(set) 
  }
  return(set)
}

make_range_svs = function(all_gr,sv_metadata_cols= c("sourceId",  "svtype", "svLen", "partner","insLen","QUAL"),
                          breakends_mean_metadata_cols=c("tumor_af","normal_af","QUAL","insLen","svLen")){
  ## last edit 2021-07-03 note replaced pairwise with faster function 

  firstStartInPair = function(gr) { start(gr) < start(partner(gr)) }
  firstInPair = function(gr) { seq_along(gr) < match(gr$partner, names(gr)) }
  #firstStartInPair=firstInPair you need this!
  
  ## Split gr into ctx and intra-chr, make pairs of intra-gr
  
  #metadata separately as for ctx both sides remain single entries
  ctx_gr = all_gr[seqnames(all_gr)!=seqnames(partner(all_gr))]
  mcols(ctx_gr) = mcols(ctx_gr)[c("bp_name",sv_metadata_cols)]
  ctx_gr$sv_name=ctx_gr$bp_name
  
  intra_gr = all_gr[seqnames(all_gr)==seqnames(partner(all_gr))]
  isFirst = firstStartInPair(intra_gr) # makes a true/false such that only select 1 partner
  #Note: if the partner range is equal, then now gets excluded -> you want to just include one of them and be done with it.
  #check if at least of of bp or partner is in the resulting array
  #otherwise, use the first in pair function
  missing_gr=intra_gr[!intra_gr$bp_name %in% c(intra_gr[isFirst]$bp_name,intra_gr[isFirst]$partner)]
  missing_isFirst = firstInPair(missing_gr)
  
  pairs =S4Vectors::Pairs(intra_gr[isFirst], partner(intra_gr)[isFirst])
  if(length(pairs)>0)  mcols(pairs)$sv_name=paste0(intra_gr[isFirst]$bp_name,"--",intra_gr[isFirst]$partner)
  
  missing_pairs =S4Vectors::Pairs(missing_gr[missing_isFirst], partner(missing_gr)[missing_isFirst])
  if(length(missing_pairs)>0) mcols(missing_pairs)$sv_name=paste0(missing_gr[missing_isFirst]$bp_name,"--",missing_gr[missing_isFirst]$partner)
  
  #append pairs and missing pairs
  pairs = c(pairs,missing_pairs)
  
  start_coord=rowMin(cbind(start(pairs@first),start(pairs@second)))
  end_coord=rowMax(cbind(end(pairs@first),end(pairs@second)))
  strand1 = pairs@first@strand %>% as.character()
  strand2 = pairs@second@strand %>% as.character()
  
  
  intra_svs = GRanges(seqnames(pairs@first),IRanges(start=start_coord,end=end_coord))
  
  metadata_pt1 = mcols(intra_gr[isFirst])[c(sv_metadata_cols,"bp_name")]
  metadata_pt2 = mcols(missing_gr[missing_isFirst])[c(sv_metadata_cols,"bp_name")]
  metadata=rbind(metadata_pt1,metadata_pt2)
  ## add bp name origin => merge with partnered bp 
  if(nrow(metadata)>0) metadata$bp_names_origin=paste0(metadata$bp_name,", ",metadata$partner)
  
  metadata$sv_name=mcols(pairs)$sv_name
  metadata$bp_name=metadata$sv_name
  
  ## in case breakends have different annotation, especially for gridss
  #breakends_mean_metadata_cols=c("tumor_af","normal_af","QUAL","insLen","svLen")
  metadata_partner_pt1 = mcols(partner(intra_gr)[isFirst])[breakends_mean_metadata_cols]
  metadata_partner_pt2 = mcols(partner(missing_gr)[missing_isFirst])[breakends_mean_metadata_cols]
  metadata_partner=rbind(metadata_partner_pt1,metadata_partner_pt2)
  
  for(col in breakends_mean_metadata_cols){
    metadata[,col] = rowMeans(cbind(metadata[,col],metadata_partner[,col]),na.rm = T)
  }
  
  
  ## name the range to both grs in the order of the new range
  ## partner also contains original bp
  metadata$partner = metadata$bp_names_origin
  mcols(intra_svs) =  metadata
  names(intra_svs) = intra_svs$bp_name
  mcols(intra_svs)$strand_origin=paste0(strand1,"/",strand2)

  ## if both partner and sv are equal in start position (and entire range) then were excluded before, now just choose 1 and add this one using its metadata   
  
  #concat ctx and intra-partnered,
  mcols(ctx_gr)$bp_names_origin = mcols(ctx_gr)$bp_name
  mcols(ctx_gr)$strand_origin = as.character(ctx_gr@strand)
  
  svs=c(intra_svs,ctx_gr)
  if(length(svs)>0)  mcols(svs)$coordinate = paste0(seqnames(svs),":",start(svs),"-",end(svs))
  
  return(svs)
}

## helper for reciprocal overlaps
overlap_fraction_gr = function(set1,set2,ignore_strand=TRUE){
  overlaps_intersect = pintersect(set1, set2,ignore.strand=ignore_strand)
  overlap_fraction = width(overlaps_intersect) / width(set1)
  return(overlap_fraction)
}

get_reciprocal_overlap_pairs = function(set1,set2,reciprocal_overlap=0.5,svtype_matching=TRUE,ignore_strand=TRUE){  
  if(is.null(names(set1))) {
    names(set1)=paste0("set1_",1:length(set1))
  }
  if(is.null(names(set2))) {
    names(set2)=paste0("set2_",1:length(set2))
  }
  # make pairs, so note that this makes the lists larger
  overlap_pairs = findOverlaps(set1,set2,type="any",minoverlap = 1,ignore.strand=ignore_strand)
  
  if(svtype_matching==TRUE){
    svtype_match = set1[overlap_pairs@from]$svtype==set2[overlap_pairs@to]$svtype
    overlap_pairs = overlap_pairs[svtype_match==TRUE]
  }
  
  #remove if not at least overlaps => reciprocal overlap
  overlap_fraction = overlap_fraction_gr(set1[overlap_pairs@from], set2[overlap_pairs@to],ignore_strand) #width of intersect/width of first
  overlap_pairs = overlap_pairs[overlap_fraction > reciprocal_overlap]
  overlap_fraction = overlap_fraction_gr(set2[overlap_pairs@to], set1[overlap_pairs@from],ignore_strand)
  overlap_pairs = overlap_pairs[overlap_fraction > reciprocal_overlap]
  
  #length of intersect divided by length of first in function 
  overlap_fraction_metadata = overlap_fraction_gr(set1[overlap_pairs@from],set2[overlap_pairs@to],ignore_strand) 
  overlap_fraction_metadata_2 = overlap_fraction_gr(set2[overlap_pairs@to],set1[overlap_pairs@from],ignore_strand)
  overlap_pairs_df = as.data.frame(overlap_pairs)
  colnames(overlap_pairs_df) = c("from","to")
    
  overlap_pairs_df$set1 = names(set1[overlap_pairs@from])
  if("svtype" %in% names(mcols(set1))) {
    overlap_pairs_df$set1_svtype = set1[overlap_pairs@from]$svtype
  }
  
  overlap_pairs_df$set2 = names(set2[overlap_pairs@to])
  if("svtype" %in% names(mcols(set2))) {
    overlap_pairs_df$set2_svtype = set2[overlap_pairs@to]$svtype
  }
  overlap_pairs_df$overlap_set1_set2 = overlap_fraction_metadata
  overlap_pairs_df$overlap_set2_set1 = overlap_fraction_metadata_2
  
  #remove overlaps with self
  overlap_pairs_df = overlap_pairs_df[overlap_pairs_df$set1!=overlap_pairs_df$set2,]
  
  return(overlap_pairs_df)
}

find_same_sv = function(set1,set2,reciprocal_overlap=0.5,svtype_matching=TRUE,ignore_strand=FALSE){  
  if(is.data.frame(set1)) { set1 =df_to_gr(set1) }
  if(is.data.frame(set2)) { set2 =df_to_gr(set2) }
  
  overlap_pairs_df = get_reciprocal_overlap_pairs(set1,set2,reciprocal_overlap=reciprocal_overlap,
                                               svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  overlap_pairs_bk = overlap_pairs_df
  #if same sets then make sure not the same overlap twice
  same_sets = FALSE
  if(all(set1$bp_name==set2$bp_name)){
    same_sets = TRUE
  }
  
  ## Merge SVs: reduce to the overlapping GRanges
  ## next map back to SV bps to group them, but make this first because sets are going to have their name reset
  # construct the dataframe from set1/2 intersection and then map those intersections to the collapsed ones to prevent multimapping rows
  
  sv_merged = GRanges()
  for(query in unique(overlap_pairs_df$from)) {
    hits = overlap_pairs_df[overlap_pairs_df$from==query,c("to")]
    if(length(hits)==0){next()}
    #sv_merged_entry = GenomicRanges::reduce(c(set1[query],set2[hits]),ignore.strand=ignore_strand)
    sv_merged_entry = GenomicRanges::intersect(set1[query],set2[hits],ignore.strand=ignore_strand)
    mcols(sv_merged_entry)$svtype= unique(set1[query]$svtype)
    sv_merged = c(sv_merged,sv_merged_entry)
    
    if(same_sets) {
      overlap_pairs_df = overlap_pairs_df[ !(overlap_pairs_df$from %in% hits & overlap_pairs_df$to==query),]
    }
  }
  
  
  if(length(sv_merged)==0) { return(data.frame())} 
  
  ## Only reduce the merged ranges that have high overlap with each other
  names(sv_merged)=paste0("merged_tmp_",1:length(sv_merged))
  
  overlap_between_merged = get_reciprocal_overlap_pairs(sv_merged,sv_merged,reciprocal_overlap=reciprocal_overlap,
                               svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  overlapping_sv_merged_names = unique(c(overlap_between_merged$set1,overlap_between_merged$set2))
  compact_sv_merged = GRanges()
  for(query in unique(overlap_between_merged$from)) {
    hits = overlap_between_merged[overlap_between_merged$from==query,c("to")]
    if(length(hits)==0){next()}
    sv_merged_entry = GenomicRanges::reduce(c(sv_merged[query],sv_merged[hits]),ignore.strand=ignore_strand)
    mcols(sv_merged_entry)$svtype= unique(sv_merged[query]$svtype)
    compact_sv_merged = c(compact_sv_merged,sv_merged_entry)
    
    overlap_between_merged = overlap_between_merged[ !(overlap_between_merged$from %in% hits & overlap_between_merged$to==query),]
  }
 
  #add the unmapped back too
  sv_merged = unique(c(compact_sv_merged,sv_merged[!names(sv_merged) %in% overlapping_sv_merged_names]))
  sv_merged = trim(sv_merged)
  names(sv_merged)=paste0("merged_",1:length(sv_merged))
  
  # Map SVs to merged
  ## require reciprocal overlap with merged range
  ## require SV type match
  
  map_sv_merged_set1_df = get_reciprocal_overlap_pairs(sv_merged,set1,reciprocal_overlap=reciprocal_overlap,
                               svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  
  
  #overlap merged to set 1 is larger than set 1 to merged, because merged ranges are larger than the original ones
  map_sv_merged_set1_df = map_sv_merged_set1_df %>% select(-to,-from) %>%
    rename(sv_merged = set1,
           sv_merged_svtype = set1_svtype,
           set1 = set2,
           set1_svtype = set2_svtype,
           overlap_merged_set1 = overlap_set1_set2,
           overlap_set1_merged = overlap_set2_set1)
  
  
  ## assign bp to sv range that is best overlapping
  map_sv_merged_set1_overlaps =  map_sv_merged_set1_df%>% group_by(set1,sv_merged) %>% summarize(overlap_mean =mean(c(overlap_merged_set1,overlap_set1_merged)))
  map_sv_merged_set1_best = map_sv_merged_set1_overlaps %>% group_by(set1) %>% summarize(max_overlap = max(as.numeric(overlap_mean)))
  map_sv_merged_set1_uq = map_sv_merged_set1_overlaps %>% merge(map_sv_merged_set1_best,by.x=c("set1","overlap_mean"), by.y=c("set1","max_overlap") )
  
  map_sv_merged_set1_df = map_sv_merged_set1_df %>% merge(map_sv_merged_set1_uq[,c("set1","sv_merged")],by=c("set1","sv_merged"))
 
  # duplicates => remove because identical
  map_sv_merged_set1_df=map_sv_merged_set1_df[!duplicated(map_sv_merged_set1_df$set1),]
  #check no duplicates map_sv_merged_set1_df[map_sv_merged_set1_df$set1 %in% map_sv_merged_set1_df[duplicated(map_sv_merged_set1_df$set1),c("set1")],]
  
  if(!same_sets){
    map_sv_merged_set2_df = get_reciprocal_overlap_pairs(sv_merged,set2,reciprocal_overlap=reciprocal_overlap,
                                                       svtype_matching=svtype_matching,ignore_strand=ignore_strand)
  
    #overlap merged to set 1 is larger than set 1 to merged, because merged ranges are larger than the original ones
    map_sv_merged_set2_df = map_sv_merged_set2_df %>% select(-to,-from) %>%
      rename(sv_merged = set1,
             sv_merged_svtype = set1_svtype,
             overlap_merged_set2 = overlap_set1_set2,
             overlap_set2_merged = overlap_set2_set1)
    
    map_sv_merged_set2_overlaps =  map_sv_merged_set2_df%>% group_by(set2,sv_merged) %>% summarize(overlap_mean =mean(c(overlap_merged_set2,overlap_set2_merged)))
    map_sv_merged_set2_best = map_sv_merged_set2_overlaps %>% group_by(set2) %>% summarize(max_overlap = max(as.numeric(overlap_mean)))
    map_sv_merged_set2_uq = map_sv_merged_set2_overlaps %>% merge(map_sv_merged_set2_best,by.x=c("set2","overlap_mean"), by.y=c("set2","max_overlap") )
    map_sv_merged_set2_df = map_sv_merged_set2_df %>% merge(map_sv_merged_set2_uq[,c("set2","sv_merged")],by=c("set2","sv_merged"))
    
    # duplicates => remove because identical
    map_sv_merged_set2_df=map_sv_merged_set2_df[!duplicated(map_sv_merged_set2_df$set2),]
    
  } else {
    map_sv_merged_set2_df = map_sv_merged_set1_df %>% rename_with(function(x){str_replace(x,"set1","set2")})
  }
  
  ## Make metadata frame
  ##map bp to overlaps and use that to build the metadata dataframe
  ## then you can merge the rows based on the reduced genomic ranges for certain properties like AF 
  # go back to the pairs that overlapped >50% amd sv type match
  
  bp_mapping = overlap_pairs_bk %>% select(-to,-from)
  
    #this second intersection is for the breakpoint mapping between the tools because does not have to be one to one 
  bp_mapping = bp_mapping %>% 
    left_join(map_sv_merged_set1_df[,c("sv_merged","sv_merged_svtype","set1","overlap_merged_set1","overlap_set1_merged")],by="set1") %>% 
    dplyr::rename(sv_merged_set1 = sv_merged)
  bp_mapping = unique(bp_mapping)
  
  bp_mapping = bp_mapping %>% 
    left_join(map_sv_merged_set2_df[,c("sv_merged","set2","overlap_merged_set2","overlap_set2_merged")],by="set2") %>%
    dplyr::rename(sv_merged_set2 = sv_merged)
  
  bp_mapping = unique(bp_mapping)
  
  ## should be identical otherwise not merged
  bp_mapping = bp_mapping %>% dplyr::mutate(sv_merged = ifelse((sv_merged_set1==sv_merged_set2),sv_merged_set1, NA))
  bp_mapping = bp_mapping %>% filter(!is.na(sv_merged)) %>% select(-sv_merged_set1,-sv_merged_set2)
  
  
  ## prevent same bp to multi ranges / select best matching group
  multimapping_merged = bp_mapping %>% group_by(sv_merged) %>%
    summarize(svs = toString(unique(sort(set1))),overlap_mean =mean(c(overlap_merged_set1,overlap_set1_merged)))
  sv_groups = multimapping_merged %>% group_by(svs) %>% summarize(max_overlap = max(overlap_mean))
  multimapping_merged_uq = multimapping_merged %>% merge(sv_groups, by.x=c("svs","overlap_mean"), by.y=c("svs","max_overlap") )
  
  #if still the same left then remove because identical
  multimapping_merged_uq = multimapping_merged_uq[!duplicated(multimapping_merged_uq$svs),]
  multimapping_merged_uq$newnames = paste0("merged_",1:nrow(multimapping_merged_uq))
  
  bp_mapping = bp_mapping %>% merge(multimapping_merged_uq[,c("newnames","sv_merged")],by=c("sv_merged")) 
  
  # add coordinates of merged
  sv_merged_coordinates = as.data.frame(sv_merged)
  sv_merged_coordinates$bp_name=row.names(sv_merged_coordinates)
  sv_merged_coordinates$sv_merged_coordinate = paste0(sv_merged_coordinates$seqnames,":",sv_merged_coordinates$start,"-",sv_merged_coordinates$end)
  if(!ignore_strand) {
    sv_merged_coordinates$sv_merged_coordinate = paste0(sv_merged_coordinates$sv_merged_coordinate,":",sv_merged_coordinates$strand)
  }
  bp_mapping= bp_mapping %>% left_join(sv_merged_coordinates[,c("bp_name","sv_merged_coordinate")],by=c("sv_merged"="bp_name"))
  
  #use newnames of selected ranges to fill out the gaps
  bp_mapping = bp_mapping %>% select(-sv_merged) %>% dplyr::rename(sv_merged = newnames)
  
  return(bp_mapping)
}

get_reciprocal_overlap_pairs_start_end = function(svs,properties,reciprocal_overlap=0,svtype_matching=F,ignore_strand=T){
  ## starting bp 
  svs_start = svs
  end(svs_start)=start(svs_start)
  
  ## ending bp 
  svs_end = svs
  start(svs_end)=end(svs_end)
  
  start_overlaps = get_reciprocal_overlap_pairs(svs_start,properties,reciprocal_overlap = reciprocal_overlap,svtype_matching = svtype_matching,ignore_strand=ignore_strand)
  if(nrow(start_overlaps)>0){
    start_overlaps$sv_breakpoint_orientation="start"
  }
  
  end_overlaps = get_reciprocal_overlap_pairs(svs_end,properties,reciprocal_overlap = reciprocal_overlap,svtype_matching = svtype_matching,ignore_strand=ignore_strand)
  if(nrow(end_overlaps)>0){
    end_overlaps$sv_breakpoint_orientation="end"
  }
  overlaps= rbind(start_overlaps,end_overlaps)
  
  return(overlaps)
}

## End of SVs as ranges 

## SV Filtering

svs_filter_chromosomes = function(svs_df,chromosomes=NULL) {
  if(is.null(chromosomes)) { 
    chromosomes= c(paste("chr",1:22,sep=""),"chrX","chrY")
  }
  chr_seqnames_df = svs_df  %>% group_by(patient_sv_merged) %>% summarize(chr_seqnames = unlist(strsplit(sv_merged_coordinate,":"))[1][])
  svs_df = svs_df %>% left_join(chr_seqnames_df,by="patient_sv_merged")
  
  svs_df = svs_df %>% filter(chr_seqnames %in% chromosomes)
  svs_df$chr_seqnames=factor(svs_df$chr_seqnames)
  return(svs_df)
}
## ENDOF filtering

## SV annotate functions

## Note: probably deprecated and replaced by long-format tables and flags 
annotate_svs_genes = function(gr,gene_properties) {
  sv_genes = findOverlaps(gr,gene_properties,ignore.strand=T)
  gr$gene_name=NA
  gr$ensembl_id=NA
  
  multiple_genes = sv_genes[duplicated(queryHits(sv_genes))]
  single_gene = sv_genes[!queryHits(sv_genes) %in% queryHits(multiple_genes)]
  
  gr[queryHits(single_gene)]$gene_name = gene_properties[subjectHits(single_gene)]$gene_name
  gr[queryHits(single_gene)]$ensembl_id = gene_properties[subjectHits(single_gene)]$gene_id
  
  for(query in unique(queryHits(multiple_genes))) {
    hits = sv_genes[queryHits(sv_genes)==query]
    gr[query]$gene_name =  toString(unique(gene_properties[subjectHits(hits)]$gene_name))
    gr[query]$ensembl_id = toString(unique(gene_properties[subjectHits(hits)]$gene_id))
  }
  
  return(gr)
}
annotate_svs_properties = function(svs,properties,cols,use_mean=F){
  properties_sv_hits = findOverlaps(properties,svs,ignore.strand=T)
  mcols(svs)[,cols]=NA
  for(query in unique(subjectHits(properties_sv_hits))) {
    hits = properties_sv_hits[subjectHits(properties_sv_hits)==query]
    hits = properties[queryHits(hits)]
    if(use_mean==F){
      mcols(svs[query])[,cols] = toString(sort(unique(mcols(hits)[,cols])))
    } else{
      mcols(svs[query])[,cols] = mean(mcols(hits)[,cols],na.rm=T)
    }
  }
  return(svs)
}
annotate_svs_features = function(gr,features_per_tx) {
  
  sv_features= findOverlaps(gr,features_per_tx,ignore.strand=T)
  gr$ft_overlap = NA
  gr$tx_id = NA
  
  multiple_features=sv_features[duplicated(queryHits(sv_features))]
  single_feature = sv_features[!queryHits(sv_features) %in% queryHits(multiple_features)]
  if(length(single_feature)>0){
    gr[queryHits(single_feature)]$ft_overlap = paste0( features_per_tx[subjectHits(single_feature)]$gene_name,":",
                                                       features_per_tx[subjectHits(single_feature)]$ft_rank, "/", features_per_tx[subjectHits(single_feature)]$ft_max)
    gr[queryHits(single_feature)]$tx_id = features_per_tx[subjectHits(single_feature)]$tx_id
  }
  
  if(length(multiple_features)>0){
    range= gr[unique(queryHits(multiple_features))]
    end(range)=start(range)
    
    sv_features_start = findOverlaps(range,features_per_tx,ignore.strand=T)
    for(query in unique(queryHits(sv_features_start))) {
      hits = sv_features_start[queryHits(sv_features_start)==query]
      
      match = features_per_tx[subjectHits(hits)]
      for(match_gene_name in unique(match$gene_name)){
        match_gene = match[match$gene_name==match_gene_name]
        if(length(unique(match_gene$tx_id))>1) { #multiple tx?? shouldt happen
          print("WARNING: multiple tx available")
          print(match_gene)
          match_gene=match_gene[match_gene$tx_id==unique(match_gene$tx_id)[1],]
        }
        if(length(match_gene)>1) { #multiple times same gene
          match_gene_ranks = paste0(min(match_gene$ft_rank),"-",max(match_gene$ft_rank))
        } else {
          match_gene_ranks = unique(match_gene$ft_rank)
        }
        
        match_feature = paste0(", start:",  match_gene_name,":",  match_gene_ranks, "/", unique(match_gene$ft_max))
        
        range[query]$ft_overlap=  paste0(range[query]$ft_overlap, match_feature)
        range[query]$tx_id = paste0(range[query]$tx_id,"start:",unique(match_gene$tx_id))
      }
    }
    gr[unique(queryHits(multiple_features))]$ft_overlap=range$ft_overlap
    gr[unique(queryHits(multiple_features))]$tx_id=range$tx_id
    
    
    
    range= gr[unique(queryHits(multiple_features))]
    start(range)=end(range)
    sv_features_end = findOverlaps(range,features_per_tx,ignore.strand=T)
    for(query in unique(queryHits(sv_features_end))) {
      hits = sv_features_end[queryHits(sv_features_end)==query]
      
      match = features_per_tx[subjectHits(hits)]
      for(match_gene_name in unique(match$gene_name)){
        match_gene = match[match$gene_name==match_gene_name]
        if(length(unique(match_gene$tx_id))>1) { #multiple tx?? shouldt happen
          print("WARNING: multiple tx available")
          print(match_gene)
          match_gene=match_gene[match_gene$tx_id==unique(match_gene$tx_id)[1],]
        }
        if(length(match_gene)>1) { #multiple times same gene
          match_gene_ranks = paste0(min(match_gene$ft_rank),"-",max(match_gene$ft_rank))
        } else {
          match_gene_ranks = unique(match_gene$ft_rank)
        }
        
        match_feature = paste0(", end:",  match_gene_name,":",  match_gene_ranks, "/", unique(match_gene$ft_max))
        
        range[query]$ft_overlap=  paste0(range[query]$ft_overlap, match_feature)
        range[query]$tx_id = paste0(range[query]$tx_id,"end:",unique(match_gene$tx_id))
      }
    }  
    
    gr[unique(queryHits(multiple_features))]$ft_overlap=range$ft_overlap
    gr[unique(queryHits(multiple_features))]$tx_id=range$tx_id
    
  }
  
  
  gr[!is.na(gr$ft_overlap)]$ft_overlap = str_replace(gr[!is.na(gr$ft_overlap)]$ft_overlap ,"NA, ", "")
  gr[!is.na(gr$tx_id)]$tx_id = str_replace(gr[!is.na(gr$tx_id)]$tx_id ,"NA, ", "")
  
  return(gr)
}

annotate_svs_exons = function(gr,gene_properties) {
  sv_genes = findOverlaps(gr,gene_properties,ignore.strand=T)
  gr$exon_gene_name=NA
  gr$exon_ensembl_id=NA
  
  multiple_genes = sv_genes[duplicated(queryHits(sv_genes))]
  single_gene = sv_genes[!queryHits(sv_genes) %in% queryHits(multiple_genes)]
  
  gr[queryHits(single_gene)]$exon_gene_name = gene_properties[subjectHits(single_gene)]$gene_name
  gr[queryHits(single_gene)]$exon_ensembl_id = gene_properties[subjectHits(single_gene)]$gene_id
  
  for(query in unique(queryHits(multiple_genes))) {
    hits = sv_genes[queryHits(sv_genes)==query]
    gr[query]$exon_gene_name =  toString(unique(gene_properties[subjectHits(hits)]$gene_name))
    gr[query]$exon_ensembl_id = toString(unique(gene_properties[subjectHits(hits)]$gene_id))
  }
  
  return(gr)
}

annotate_svs_onco_tsg = function(svs_anno_df,ensembl_id_colname="ensembl_id",ensembl_gene_map,oncogenes,tsg){
  #patient, sv_name, separate on annotated genes => rows, then merge gene name with dataframe / ensembl ID needed
  
  select_cols = c("bp_name",ensembl_id_colname)
  svs_ensembl_id_long = separate_rows(svs_anno_df[,select_cols], all_of(ensembl_id_colname), sep=", ")
  svs_ensembl_id_long[,c("ensembl_id")] = lapply(svs_ensembl_id_long[,ensembl_id_colname],remove_version_from_id)
  
  svs_ensembl_id_long = svs_ensembl_id_long %>% left_join(ensembl_gene_map,by=c("ensembl_id"))
  svs_ensembl_id_long = svs_ensembl_id_long %>% dplyr::rename(ensembl_id_version = gene_id)
  
  svs_ensembl_id_long = svs_ensembl_id_long %>% 
    dplyr::mutate( oncogene= (ensembl_id %in% oncogenes |
                                gene_name %in% oncogenes)) %>%
    dplyr::mutate( tsg = (ensembl_id %in% tsg |
                            gene_name %in% tsg))
  
  
  svs_anno_df = svs_anno_df %>% 
    dplyr::mutate(has_oncogene = bp_name %in% svs_ensembl_id_long[svs_ensembl_id_long$oncogene,]$bp_name) %>%
    dplyr::mutate(has_tsg = bp_name %in% svs_ensembl_id_long[svs_ensembl_id_long$tsg,]$bp_name) 
  return(svs_anno_df)
}

## End of SV annotate functions

## Transcript selection 

get_tx_by_genes = function(gtf,ensembl_id_lst) {
  transcripts = gtf[gtf$type=="transcript"& gtf$gene_id %in% ensembl_id_lst]
  metadata_tx = as.data.frame(mcols(transcripts))
  metadata_tx$mane_select=F
  metadata_tx[metadata_tx$transcript_id %in% tx_mane_select,c("mane_select")]=TRUE
  
  protein_max_exon = as.data.frame(gtf[!is.na(gtf$protein_id)&gtf$gene_id %in% ensembl_id_lst]) %>% group_by(protein_id) %>%
    summarize(max_exon = max(as.numeric(exon_number),na.rm=T))
  protein_max_cds = as.data.frame(gtf[!is.na(gtf$protein_id)&gtf$gene_id %in% ensembl_id_lst & gtf$type=="CDS"]) %>% group_by(protein_id) %>%
    summarize(max_cds_width = sum(width,na.rm=T))
  
  metadata_tx = metadata_tx %>% left_join(protein_max_exon,by="protein_id")
  metadata_tx = metadata_tx %>% left_join(protein_max_cds,by="protein_id")
  
  mcols(transcripts)=metadata_tx
  return(transcripts)
}

filter_tx_canonical = function(tx_table) {
  
  if(length(tx_table[tx_table$mane_select])>0) {
    tx_table = tx_table[tx_table$mane_select]
  }
  
  if(length(tx_table[grepl("basic|CCDS|apris",tx_table$tag)])>0) {
    tx_table = tx_table[grepl("basic|CCDS",tx_table$tag)]
  }
  
  if(length(tx_table[tx_table$transcript_type=="protein_coding",])>0) {
    tx_table = tx_table[tx_table$transcript_type=="protein_coding"]
  }
  
  if(length(tx_table)==1) return(tx_table)
  
  max_cds_width = as.numeric(max(tx_table$max_cds_width,na.rm=T))
  if(length(tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width > max_cds_width*0.9])>0) {
    tx_table = tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width > max_cds_width*0.9]
  }
  
  if(length(tx_table)==1) return(tx_table)
  
  max_exons = as.numeric(max(tx_table$max_exon,na.rm=T))*0.8
  if(length(tx_table[!is.na(tx_table$max_exons) & tx_table$max_exons > max_exons])>0) {
    tx_table = tx_table[!is.na(tx_table$max_exons) & tx_table$max_exons > max_exons]
  }
  
  if(length(tx_table)==1) return(tx_table)
  
  min_tsl = min(tx_table$transcript_support_level,na.rm=T)
  if(length(tx_table[!is.na(tx_table$transcript_support_level) & tx_table$transcript_support_level == min_tsl])>0) {
    tx_table = tx_table[!is.na(tx_table$transcript_support_level) & tx_table$transcript_support_level == min_tsl]
  }
  
  
  if(length(tx_table)==1) return(tx_table)
  
  tx_table$tx_width = width(tx_table@ranges)
  max_width = max(tx_table$tx_width,na.rm=T)
  if(length(tx_table[tx_table$tx_width == max_width])>0) {
    tx_table = tx_table[tx_table$tx_width == max_width]
  }
  mcols(tx_table) = mcols(tx_table)[names(mcols(tx_table)) != "tx_width"]
  
  
  if(length(tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width == max_cds_width])>0) {
    tx_table = tx_table[!is.na(tx_table$max_cds_width) & tx_table$max_cds_width == max_cds_width]
  }
  
  return(tx_table)
}

get_features_per_tx = function(features, transcript_lst){
  
  features_per_tx=GRanges()
  #use names because sometimes empty
  for(tx_id in transcript_lst) {
    tx_features = features[[tx_id]]
    if(length(tx_features)==0){next()}
    tx_features$tx_id =tx_id
    
    tx_features$ft_rank = rank(tx_features)
    tx_features$ft_max = max(tx_features$ft_rank,na.rm=T)
    if(any(strand(tx_features)=="-")) {
      tx_features$ft_rank = (tx_features$ft_max-tx_features$ft_rank)+1
    }
  
    features_per_tx = c(features_per_tx,tx_features)
  }
  return(features_per_tx)
}

## SV burden:
cnt_svs = function(svs_df, attr_name=NULL, group_cols=c("patient_id")) {
  cnt_table = svs_df %>% group_by(across(all_of(group_cols))) %>% summarise(svs_cnt = n(),.groups="keep") %>% ungroup()
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c(group_cols,attr_name)
  }
  return(cnt_table)
}


## Loading can be standardized across utils scripts
read_svs_df = function(svs_path,patient_id=NULL) {
  svs_df = read.table(svs_path,header=T,sep="\t",stringsAsFactors=F) 
  
  sv_columns = names(svs_df)
  
  if(!"patient_sv_name" %in% sv_columns) {
    if(!"sv_name" %in% sv_columns) {
      print("ERROR: missing patient_sv_name and sv_name")
      quit()
    }
    if(!"patient_id" %in% sv_columns) {
      if(is.null(patient_id)) {
        print("WARNING: missing patient_id, proceeding with sv_name instead of patient_sv_name")
        svs_df$patient_sv_name = svs_df$sv_name
      } else {
        svs_df$patient_id = patient_id
        svs_df$patient_sv_name = paste0(svs_df$patient_id,"_",svs_df$sv_name)
      }
    } 
  }
  
  #for reporting
  svs_df=svs_df %>% filter(grepl("chr",coordinate))
  svs_df$from_coordinate = svs_df$coordinate
  return(svs_df)
}

get_svs_gr = function(svs_df) {
  svs_gr = GRanges(svs_df$coordinate)
  svs_gr$patient_sv_name=svs_df$patient_sv_name
  svs_gr$svtype=svs_df$svtype
  names(svs_gr)=svs_gr$patient_sv_name
  return(svs_gr)
}

## used in sv burden
load_patient_svs = function(sv_patient_path,patient) {
  sv_patient = read.table(sv_patient_path,header = T,sep="\t",stringsAsFactors = F)
  sv_patient$patient_id=patient$patient_id
  
  sv_patient$patient_sv_merged = paste0(sv_patient$patient_id,"_",sv_patient$sv_merged)
  sv_patient$patient_sv_name = paste0(sv_patient$patient_id,"_",sv_patient$sv_name)
  sv_patient$partner_sv_name = paste0(sv_patient$patient_id,"_",sv_patient$partner)
  
  
  #do NOT filter 
  #sv_patient = svs_filter_chromosomes(sv_patient)
  sv_patient = annotate_sv_af_class(sv_patient)
  
  return(sv_patient)
}


get_multi_tool_support = function(svs_df){
  multi_tool_support = svs_df %>% filter(grepl("merged",patient_sv_merged)) %>% group_by(patient_sv_merged,svtype) %>% summarize(tool_cnt=length(unique(tool)))
  multi_tool_support = filter(multi_tool_support,tool_cnt>1)
  
  # For CTX: require partners to be supported too 
  if(is.null(svs_df$partner_sv_name) & !is.null(svs_df$patient_id)){
    svs_df$partner_sv_name = paste0(svs_df$patient_id,"_",svs_df$partner)
  }
  
  multi_tool_ctx = filter(svs_df, patient_sv_merged %in% filter(multi_tool_support,tool_cnt>1&svtype=="CTX")$patient_sv_merged)
  ctx_partners = svs_df %>% filter(patient_sv_name %in% multi_tool_ctx$partner_sv_name & patient_sv_merged %in% multi_tool_ctx$patient_sv_merged)
  
  #svs filtered out
  #svs_df %>% filter((patient_sv_name %in% multi_tool_ctx$partner_sv_name & !patient_sv_merged %in% multi_tool_ctx$patient_sv_merged) | 
  #                    (patient_sv_merged %in% multi_tool_ctx$patient_sv_merged & !partner_sv_name %in% multi_tool_ctx$patient_sv_name))
  
  multi_tool_support=multi_tool_support %>% filter(svtype!="CTX" |
                                                     patient_sv_merged %in%  ctx_partners$patient_sv_merged)
  #usage: filter(patient_sv_merged %in% multi_tool_support$patient_sv_merged)
  
  return(multi_tool_support)
}


## SV overlaps to genes

make_sv_bp_gene_partnered = function(sv_bp_gene_overlaps, gene_overlap_cols=c("gene_type","gene_name","gene_id"),svs_df_cols = c("patient_sv_name","patient_sv_merged","partner_sv_name","svtype")) {
  
  svs_unique_cols=c(svs_df_cols,gene_overlap_cols)
  svs_unique_cols=svs_unique_cols[svs_unique_cols %in% names(sv_bp_gene_overlaps)]
  
  svs_to_genes = sv_bp_gene_overlaps[,svs_unique_cols] %>% unique()
  
  
  ## For CTX only => rest will not have partner sv name matching like this, 
  # but do not use     filter(patient_sv_name %in% svs_to_genes$partner_sv_name) because then no longer svs with 1 edge 
  #if you dont remove sv bp orientation you will get inflated counts
  
  sv_bp_gene_partnered_ctx = svs_to_genes %>% filter(svtype=="CTX") %>%
    left_join(svs_to_genes %>% select(-svtype),
              by=c("partner_sv_name"="patient_sv_name","patient_sv_name"="partner_sv_name"))
  
  sv_bp_gene_partnered_ctx = sv_bp_gene_partnered_ctx %>% 
    dplyr::rename_with(.cols = ends_with(".y"), function(x){paste0("partner_",substr(x,0,stop = (str_length(x)-2)))}) %>% 
    dplyr::rename_with(.cols = ends_with(".x"), function(x){substr(x,0,stop = (str_length(x)-2))})
  
  sv_bp_gene_partnered_ctx = sv_bp_gene_partnered_ctx %>% dplyr::rename(partner_sv_merged = partner_patient_sv_merged)
  
  #with sv bp orientation for intra-svs
  svs_to_genes_orientation = sv_bp_gene_overlaps[,c(svs_unique_cols,"sv_breakpoint_orientation")] %>% unique()
  
  #in case only 'end' coordinate, include it in the first datafame and will not have a partner in the 2nd
  sv_bp_gene_partnered_intra = svs_to_genes_orientation %>% filter(svtype!="CTX") %>%
    filter(sv_breakpoint_orientation=="start" | !patient_sv_name %in% filter(svs_to_genes_orientation,sv_breakpoint_orientation=="start")$patient_sv_name ) %>% 
    select(-sv_breakpoint_orientation) %>% unique() 
  sv_bp_gene_partnered_intra_2 = svs_to_genes_orientation %>% filter(svtype!="CTX") %>%
    filter(sv_breakpoint_orientation=="end" & patient_sv_name %in% filter(svs_to_genes_orientation,sv_breakpoint_orientation=="start")$patient_sv_name ) %>% 
    select(patient_sv_name,all_of(gene_overlap_cols)) %>% unique() %>%
    dplyr::rename_with(.cols = all_of(gene_overlap_cols), function(x){paste0("partner_",x)}) 
  
  sv_bp_gene_partnered_intra = sv_bp_gene_partnered_intra %>% left_join(sv_bp_gene_partnered_intra_2)
  
  sv_bp_gene_partnered_intra$partner_sv_merged = sv_bp_gene_partnered_intra$patient_sv_merged
  
  sv_bp_gene_partnered_cols = c(names(sv_bp_gene_partnered_ctx),names(sv_bp_gene_partnered_intra)) %>% unique()
  
  sv_bp_gene_partnered = rbind(sv_bp_gene_partnered_ctx[,sv_bp_gene_partnered_cols],
                               sv_bp_gene_partnered_intra[,sv_bp_gene_partnered_cols])
  
  
  #check if all rows are there
  #nrow(sv_bp_gene_overlaps %>% filter(!patient_sv_name %in% sv_bp_gene_partnered$patient_sv_name))!=0
  
  return(sv_bp_gene_partnered)
} 


make_sv_bp_gene_partnered = function(sv_bp_gene_overlaps, 
                                            gene_overlap_cols=c("gene_type","gene_name","gene_id"),
                                            svs_df_cols = NULL,
                                            patient_sv_col="patient_sv_name",partner_sv_col="partner_sv_name") {
#patient_sv_col="patient_sv_merged";partner_sv_col="partner_sv_merged"
  if(is.null(svs_df_cols)) {
    svs_df_cols = c(patient_sv_col,partner_sv_col,"patient_sv_merged","svtype") %>% unique()
  }
  svs_unique_cols=c(svs_df_cols,gene_overlap_cols)
  svs_unique_cols=svs_unique_cols[svs_unique_cols %in% names(sv_bp_gene_overlaps)]
  
  svs_to_genes = sv_bp_gene_overlaps[,svs_unique_cols] %>% unique()
  
  
  ## For CTX only => rest will not have partner sv name matching like this, 
  # but do not use     filter(patient_sv_name %in% svs_to_genes$partner_sv_name) because then no longer svs with 1 edge 
  #if you dont remove sv bp orientation you will get inflated counts
  
  sv_bp_gene_partnered_ctx = svs_to_genes %>% filter(svtype=="CTX") %>%
    left_join(svs_to_genes %>% select(-svtype)%>% dplyr::rename(!!sym(patient_sv_col) := !!sym(partner_sv_col),!!sym(partner_sv_col) := !!sym(patient_sv_col)),
              by=c(partner_sv_col,patient_sv_col))
  #by=c(!!sym(partner_sv_col) := !!sym(patient_sv_col), !!sym(patient_sv_col) := !!sym(partner_sv_col)))
  
  sv_bp_gene_partnered_ctx = sv_bp_gene_partnered_ctx %>% 
    dplyr::rename_with(.cols = ends_with(".y"), function(x){paste0("partner_",substr(x,0,stop = (str_length(x)-2)))}) %>% 
    dplyr::rename_with(.cols = ends_with(".x"), function(x){substr(x,0,stop = (str_length(x)-2))})
  
  if(partner_sv_col!="partner_sv_merged" & "partner_patient_sv_merged" %in% names(sv_bp_gene_partnered_ctx)) {
  sv_bp_gene_partnered_ctx = sv_bp_gene_partnered_ctx %>% dplyr::rename(!!sym(partner_sv_col) := partner_patient_sv_merged)
  }
  
  #with sv bp orientation for intra-svs
  svs_to_genes_orientation = sv_bp_gene_overlaps[,c(svs_unique_cols,"sv_breakpoint_orientation")] %>% unique()
  
  #in case only 'end' coordinate, include it in the first datafame and will not have a partner in the 2nd
  sv_bp_gene_partnered_intra = svs_to_genes_orientation %>% filter(svtype!="CTX") %>%
    filter(sv_breakpoint_orientation=="start" | ! (!!sym(patient_sv_col)) %in% filter(svs_to_genes_orientation,sv_breakpoint_orientation=="start")[,patient_sv_col] ) %>% 
    select(-sv_breakpoint_orientation) %>% unique() 
  sv_bp_gene_partnered_intra_2 = svs_to_genes_orientation %>% filter(svtype!="CTX") %>%
    filter(sv_breakpoint_orientation=="end" & !!sym(patient_sv_col) %in% filter(svs_to_genes_orientation,sv_breakpoint_orientation=="start")[,patient_sv_col] ) %>% 
    select(all_of(c(patient_sv_col,gene_overlap_cols))) %>% unique() %>%
    dplyr::rename_with(.cols = all_of(gene_overlap_cols), function(x){paste0("partner_",x)}) 
  
  sv_bp_gene_partnered_intra = sv_bp_gene_partnered_intra %>% left_join(sv_bp_gene_partnered_intra_2)
  
  
  sv_bp_gene_partnered_intra[,partner_sv_col] = sv_bp_gene_partnered_intra[,patient_sv_col]
  
  sv_bp_gene_partnered_cols = c(names(sv_bp_gene_partnered_ctx),names(sv_bp_gene_partnered_intra)) %>% unique()
  
  sv_bp_gene_partnered = rbind(sv_bp_gene_partnered_ctx[,sv_bp_gene_partnered_cols],
                               sv_bp_gene_partnered_intra[,sv_bp_gene_partnered_cols])
  
  
  #check if all rows are there
  #nrow(sv_bp_gene_overlaps %>% filter(!patient_sv_merged %in% sv_bp_gene_partnered$patient_sv_merged))!=0
  
  return(sv_bp_gene_partnered)
}

## generic annotation used in sv burden but also recurrence analyses

make_svlen_bin=function(svs_df){
  svs_df= svs_df %>% mutate(svlen_bin = ifelse(is.na(svLen),"",
                                               ifelse(svLen<1000,"<1kb",
                                                      ifelse(svLen<10000,"<10kb",
                                                             ifelse(svLen<500000,"<500kb",
                                                                    ifelse(svLen<5000000,"<5Mb",">5Mb"))))))
  return(svs_df)
}


make_tumor_af_bin=function(svs_df){
  svs_df= svs_df %>% mutate(tumor_af_bin = ifelse(is.na(tumor_af),"",
                                                  ifelse(tumor_af<0.1,"<0.1taf",
                                                         ifelse(tumor_af<0.2,"0.1<taf<0.2",">0.2taf"))))
  return(svs_df)
}

## annotation in sv overlaps context eg with genes exons or repeats

annotate_sv_affects_exon = function(df) {
  df = df %>% mutate(flag_affects_exon = (svtype=="CTX" & flag_bp_in_gene==T) | flag_overlap_exon==T | flag_bp_in_exon==T)
  return(df)
}

annotate_sv_repeat_mediated = function(df) {
  df = df %>% mutate(flag_repeat_mediated = (flag_repeat_both_bp==T | flag_segdup_both_bp==T))
  return(df)
}

annotate_sv_gene_position = function(gene_centric_df) {
  ## todo make more generic / set1_set2 etc
  gene_centric_df = gene_centric_df %>% mutate(sv_position = 
                                                 ifelse(overlap_gene_sv_frac>0.99,"spanning",
                                                        ifelse(overlap_sv_gene_frac>0.99,"inside",
                                                               "partial_overlap")))
  
  return(gene_centric_df)
}
##


flag_sv_overlap = function(svs_df,overlaps,sv_key_col="patient_sv_name",overlaps_key_col="set1",overlaps_colname="overlap") {
  svs_df[,overlaps_colname]=F
  svs_df[svs_df[,sv_key_col] %in% overlaps[,overlaps_key_col], overlaps_colname]=T
  return(svs_df)
}

#function to receive overlap attributes, e.g. identifier and overlap%
annotate_sv_overlap_attr = function(svs_df,overlaps,attr_col=NULL,attr_max_cnt=5,report_attr_cnt=F,report_overlap_frac=T,sv_key_col="patient_sv_name",overlaps_key_col="set1",overlaps_colname="overlap") {
  ## todo add overlap frac, but can be multiple entries (.e.g genes) and should report for every single one which requires merging over the svs first 
  
  if(is.null(attr_col)) {
    attr_col="attr_col"
  }
  
  if(report_overlap_frac){
    summary_overlaps = overlaps %>% group_by(across(all_of(overlaps_key_col))) %>% 
      summarize(
        attr_cnt=length(unique(.data[[attr_col]])),
        attr_summary = ifelse(attr_cnt <= attr_max_cnt, 
                              toString(unique(sort(
                                paste0(.data[[attr_col]]," (",round(overlap_set1_set2,2),", ",round(overlap_set2_set1,2),")")
                              ))),
                              paste0(attr_cnt,">",attr_max_cnt," ",attr_col) ) )
    
  } else {
    summary_overlaps = overlaps %>% group_by(across(all_of(overlaps_key_col))) %>% 
      summarize(
        attr_cnt=length(unique(.data[[attr_col]])),
        attr_summary = ifelse(attr_cnt <= attr_max_cnt, 
                              toString(unique(sort(.data[[attr_col]]))),
                              paste0(attr_cnt,">",attr_max_cnt," ",attr_col) ) )
  }
  
  summary_overlaps = summary_overlaps %>% 
    dplyr::rename(!!overlaps_colname := attr_summary) %>% ungroup() %>% as.data.frame()
  
  if(report_attr_cnt==F) {
    summary_overlaps = summary_overlaps %>% select(-attr_cnt)
  } else {
    overlaps_colname_cnt=paste0(overlaps_colname,"_cnt")
    summary_overlaps = summary_overlaps %>% dplyr::rename(!!overlaps_colname_cnt := attr_cnt)
  }
  
  svs_df = svs_df %>% merge(summary_overlaps, all.x=T, by.x =sv_key_col, by.y= overlaps_key_col)
  
  return(svs_df)
}

annotate_sv_cytoband = function(svs_df,chromosome_bands) {
  chromosome_bands = chromosome_bands[as.character(chromosome_bands@seqnames) %in% autosomes]
  
  if("cytoband" %in% names(svs_df)) {
    print("Error: contains cytoband column, remove first.")
    return(svs_df)
  }
  
  svs = GRanges(svs_df)
  names(svs) = svs$patient_sv_name
  
  svs_cytoband_overlaps = get_reciprocal_overlap_pairs_start_end(svs,chromosome_bands,reciprocal_overlap = 0,svtype_matching = F)
  
  svs_cytoband_overlaps = svs_cytoband_overlaps %>% select(set1,sv_breakpoint_orientation,set2) %>% pivot_wider(names_from="sv_breakpoint_orientation",values_from="set2")
  svs_cytoband_overlaps[is.na(svs_cytoband_overlaps)]=""
  svs_cytoband_overlaps =  svs_cytoband_overlaps %>% dplyr::rename(cytoband_start=start,cytoband_end=end,patient_sv_name=set1)
  svs_cytoband_overlaps = svs_cytoband_overlaps %>% group_by(patient_sv_name) %>% 
    summarize(cytoband = ifelse(cytoband_start!=cytoband_end, 
                                paste0(cytoband_start,"-",cytoband_end),
                                cytoband_start))
  
  svs_df = svs_df %>% left_join(svs_cytoband_overlaps[,c("patient_sv_name","cytoband")], by=c("partner_sv_name"="patient_sv_name")) %>%
    dplyr::rename(partner_cytoband=cytoband) %>% left_join(svs_cytoband_overlaps[,c("patient_sv_name","cytoband")], by=c("patient_sv_name")) 
  
  
  return(svs_df)
}
annotate_sv_chr_arm = function(svs_df,chr_arms) {
  
  if("chr_arm" %in% names(svs_df)) {
    print("Error: contains chr_arm column, remove first.")
    return(svs_df)
  }
  
  
  svs = GRanges(svs_df)
  names(svs) = svs$patient_sv_name
  
  svs_chr_arm_overlaps = get_reciprocal_overlap_pairs_start_end(svs,chr_arms,reciprocal_overlap = 0,svtype_matching = F)
 
  svs_chr_arm_overlaps = svs_chr_arm_overlaps %>% select(set1,sv_breakpoint_orientation,set2) %>%  pivot_wider(names_from="sv_breakpoint_orientation",values_from="set2")
  svs_chr_arm_overlaps[is.na(svs_chr_arm_overlaps)]=""
  svs_chr_arm_overlaps =  svs_chr_arm_overlaps %>% dplyr::rename(chr_arm_start=start,chr_arm_end=end,patient_sv_name=set1)
  svs_chr_arm_overlaps = svs_chr_arm_overlaps %>% group_by(patient_sv_name) %>% 
    summarize(chr_arm = ifelse(chr_arm_start!=chr_arm_end, 
                                paste0(chr_arm_start,"-",chr_arm_end),
                                chr_arm_start))
  
  
  svs_df = svs_df %>% left_join(svs_chr_arm_overlaps[,c("patient_sv_name","chr_arm")], by=c("partner_sv_name"="patient_sv_name")) %>%
    dplyr::rename(partner_chr_arm=chr_arm) %>% left_join(svs_chr_arm_overlaps[,c("patient_sv_name","chr_arm")], by=c("patient_sv_name")) 
  
  
  return(svs_df)
}



get_sv_gene_flank_overlaps = function(svs,genes_1mb,svs_df,gene_properties_df){
  svs_df_overlap_cols = svs_df_overlap_cols[svs_df_overlap_cols %in% names(svs_df)]
  svs_df_anno_cols = svs_df_anno_cols[svs_df_anno_cols %in% names(svs_df)]
  gene_properties_df_cols = gene_properties_df_cols[gene_properties_df_cols %in% names(gene_properties_df)]
  gene_properties_df_flags = gene_properties_df_flags[gene_properties_df_flags %in% names(gene_properties_df)]
  
  
  flank_gene_overlaps= get_reciprocal_overlap_pairs_start_end(svs,genes_1mb,reciprocal_overlap = 0,svtype_matching = FALSE)
  
  join_gene_properties_df = gene_properties_df[,c("strand",gene_properties_df_cols,gene_properties_df_flags,"cytoband","start","end")] %>% dplyr::rename(gene_start=start,gene_end=end)
  
  join_svs_df = svs_df[,c(svs_df_overlap_cols,svs_df_anno_cols,"start","end")] %>% dplyr::rename(sv_start=start,sv_end=end)
  join_svs_df = unique(join_svs_df)
  
  flank_gene_overlaps = flank_gene_overlaps %>% dplyr::rename(patient_sv_name=set1, gene_id=set2, svtype=set1_svtype) %>% 
    left_join(join_gene_properties_df)  %>%
    left_join(join_svs_df) 
  
  if(exists("annotate_gene_group")){
  flank_gene_overlaps = flank_gene_overlaps %>% annotate_gene_group(cancer_genes)
  }
  #sv start always << end 
  ## upstream if start & end << gene coordinate and gene is on + strand 
  
  
  #for full sv/gene so spannng becomes 0 flank_gene_overlaps$distance = distance(GRanges(flank_gene_overlaps$to_coordinate),GRanges(flank_gene_overlaps$from_coordinate))
  
  #check 
  flank_gene_overlaps = flank_gene_overlaps %>% rowwise() %>% mutate(distance = ifelse(sv_breakpoint_orientation=="start",
                                                                                       min(c(abs(sv_start-gene_start),abs(sv_start-gene_end))),
                                                                                       min(c(abs(sv_end-gene_start),abs(sv_end-gene_end))))
  ) %>% as.data.frame() 
  flank_gene_overlaps$distance_kbp = flank_gene_overlaps$distance/1e3
  return(flank_gene_overlaps) 
}

get_sv_cluster_mapping = function(graph) {
  independent_sets= igraph::clusters(graph,mode="strong")
  sv_cluster_mapping=independent_sets$membership %>% as.data.frame()
  colnames(sv_cluster_mapping)=c("cluster")
  sv_cluster_mapping$patient_sv_merged = rownames(sv_cluster_mapping)
  
  return(sv_cluster_mapping)
}


get_orphan_partner_svs_df = function(orphan_source_svs_df) {
  #orphan_source_svs_df %>% select(patient_sv_merged,patient_sv_name,coordinate,partner_sv_name,partner_coordinate)
  #make own metadata df => switch partner and patient for processing
  orphan_partner_svs = GRanges(orphan_source_svs_df$partner_coordinate) 
  names(orphan_partner_svs)=orphan_source_svs_df$partner_sv_name
  orphan_partner_svs$svtype="CTX"
  orphan_partner_svs$patient_sv_name=names(orphan_partner_svs)
  orphan_partner_svs$partner_sv_merged=orphan_source_svs_df$patient_sv_merged
  orphan_partner_svs$partner_sv_name=orphan_source_svs_df$patient_sv_name
  orphan_partner_svs$patient_label=orphan_source_svs_df$patient_label
  
  orphan_partner_svs_df = orphan_partner_svs %>% data.frame() %>% get_gr_coordinate()
  orphan_partner_svs_df$coordinate = paste0(orphan_partner_svs_df$coordinate,":",orphan_partner_svs_df$strand)
  return(orphan_partner_svs_df)
}
get_orphan_partner_merged = function(orphan_merged_svs, unfiltered_svs_df) {
  library(igraph)
  library(openssl)
  
  
  orphan_source_svs_df = unfiltered_svs_df %>% filter(patient_sv_merged %in% filter(orphan_merged_svs)$patient_sv_merged) %>% unique()
  orphan_partner_svs_df = get_orphan_partner_svs_df(orphan_source_svs_df)
  orphan_partner_svs = GRanges(orphan_partner_svs_df)
  names(orphan_partner_svs) = orphan_partner_svs$patient_sv_name
  
  orphan_distance=1e6
  
  partner_overlaps = get_reciprocal_overlap_pairs(resize_gr_distance(orphan_partner_svs,orphan_distance),resize_gr_distance(orphan_partner_svs,orphan_distance),ignore_strand = F)
  
  orphan_partner_svs_cols = names(orphan_partner_svs_df)
  orphan_partner_svs_cols = orphan_partner_svs_cols[!orphan_partner_svs_cols %in% c("svtype")]
  orphan_partner_svs_cols_display = orphan_partner_svs_cols[!orphan_partner_svs_cols %in% c("patient_sv_name")]
  
  partner_overlaps = partner_overlaps %>% 
    left_join(orphan_partner_svs_df[,orphan_partner_svs_cols],by=c("set1"="patient_sv_name")) %>% dplyr::rename_with(.cols=orphan_partner_svs_cols_display,.fn=function(x){paste0("set1_",x)}) %>%
    left_join(orphan_partner_svs_df[,orphan_partner_svs_cols],by=c("set2"="patient_sv_name")) %>% dplyr::rename_with(.cols=orphan_partner_svs_cols_display,.fn=function(x){paste0("set2_",x)}) 
  
  #same patient, same partner sv merged, but not same sv
  partner_overlaps = partner_overlaps %>% 
    filter(set1_patient_label==set2_patient_label) %>%
    filter(set1!=set2) %>%
    filter(set1_partner_sv_merged==set2_partner_sv_merged) %>% unique()
  
  partner_overlaps$distance = GenomicRanges::distance(GRanges(partner_overlaps$set1_coordinate),GRanges(partner_overlaps$set2_coordinate),ignore.strand=T)
  
  graph = graph_from_data_frame(partner_overlaps[partner_overlaps$distance<1e3,c("set1","set2")],directed=F)
  
  cliques=maximal.cliques(graph)
  orphan_cluster_mapping = get_sv_cluster_mapping(graph) %>% dplyr::rename(patient_sv_name=patient_sv_merged)
  orphan_partner_sv_clusters = orphan_partner_svs_df %>% left_join(orphan_cluster_mapping) %>%  group_by(partner_sv_merged,cluster) %>% summarize(svs=toString(sort(unique(patient_sv_name))))
  
  ## make merged id hash to fill slot of unfiltered svs and get them through merging
  orphan_partner_sv_clusters = orphan_partner_svs_df %>% left_join(orphan_cluster_mapping) %>% 
    filter(!is.na(cluster)) %>% 
    group_by(partner_sv_merged,cluster) %>% summarize(svs_lst=toString(sort(unique(patient_sv_name))), orphan_merged_id=md5(svs_lst),
                                                      seqnames=unique(seqnames),start=min(start),end=max(end),strand=unique(strand)) %>% ungroup() %>% as.data.frame()
  #remove the ones that are NA cluster otherwise all not overlapping will be 1 group
  orphan_partner_sv_clusters = orphan_partner_sv_clusters %>%  filter(!is.na(cluster)) 
  
  if(nrow(orphan_partner_sv_clusters)>0) {
    
    orphan_partner_sv_clusters = orphan_partner_sv_clusters %>% get_gr_coordinate("partner_sv_merged_coordinate")
    orphan_partner_sv_clusters$partner_sv_merged_coordinate = paste0(orphan_partner_sv_clusters$partner_sv_merged_coordinate,":",orphan_partner_sv_clusters$strand)
    
    #ifelse to prevent M..._orphan_merged_NA
    orphan_partner_sv_clusters = orphan_partner_sv_clusters %>% select(partner_sv_merged,orphan_merged_id,partner_sv_merged_coordinate,)
    orphan_source_svs_df = orphan_source_svs_df %>% left_join(orphan_partner_sv_clusters,by=c("patient_sv_merged"="partner_sv_merged")) %>%
      dplyr::mutate(partner_sv_merged=ifelse(is.na(orphan_merged_id),NA,paste0(patient_label,"_orphan_merged_",orphan_merged_id))) %>% 
      select(-orphan_merged_id) %>%
      as.data.frame() %>% unique() 
  } else {
    orphan_source_svs_df$partner_sv_merged=NA
    orphan_source_svs_df$partner_sv_merged_coordinate=NA
  }
  
  
  #orphan_source_svs_df %>%
  #  select(patient_label,patient_sv_merged,patient_sv_name,coordinate,partner_sv_name,partner_coordinate,partner_sv_merged) %>% View()
  
  return(orphan_source_svs_df)
}

resolve_orphan_svs = function(merged_svs,unfiltered_svs_df,as_global=T) {
  #one end merged and check other end
  #the partner bps not merged => separate function, is that really needed? Maybe not but easier because only filtered
  orphan_merged_svs = merged_svs %>% filter(svtype=="CTX") %>% filter(is.na(partner_sv_merged_coordinate)) %>% unique()
  if(orphan_merged_svs %>% nrow()==0) { return(T)}
  orphan_source_svs_df = get_orphan_partner_merged(orphan_merged_svs,unfiltered_svs_df)
  
  #these svs have source in merged now better annotated
  orphan_merged_svs = orphan_merged_svs %>% select(-partner_sv_merged,-partner_sv_merged_coordinate) %>% left_join(orphan_source_svs_df[,c("patient_sv_merged","partner_sv_merged","partner_sv_merged_coordinate")] %>% unique())
  #remove then add back, but only the ones that do not have NA of course
  remain_orphan_svs = unfiltered_svs_df %>% filter(patient_sv_merged %in% filter(orphan_merged_svs,is.na(partner_sv_merged))$patient_sv_merged)
  orphan_merged_svs = orphan_merged_svs %>% filter(!is.na(partner_sv_merged))
  
  if(remain_orphan_svs %>% nrow() > 0) {
    print("These merged svs had non matching partners")
    remain_orphan_svs %>% unique()
  }
  merged_svs = merged_svs %>% filter(ifelse(svtype=="CTX",!is.na(sv_merged_coordinate)&!is.na(partner_sv_merged_coordinate),TRUE))
  merged_svs = rbind(merged_svs,orphan_merged_svs)
  
  #should also add the partner orphan as 'merged' such that both sides are present in the merged sv df
  #this should exist: merged_svs %>% filter(patient_sv_merged %in% orphan_merged_svs$partner_sv_merged)
  
  orphan_partner_svs_df = get_orphan_partner_svs_df(orphan_source_svs_df) %>% 
    left_join(orphan_source_svs_df[,c("partner_sv_name","partner_sv_merged","partner_sv_merged_coordinate")] %>%
                dplyr::rename(patient_sv_merged=partner_sv_merged,patient_sv_merged_coordinate=partner_sv_merged_coordinate),
              by=c("patient_sv_name"="partner_sv_name")) 
  
  #the partner ends with their new merged ids.
  orphan_partner_svs_df = orphan_partner_svs_df %>% filter(patient_sv_merged %in% orphan_merged_svs$partner_sv_merged) %>% 
    select(c("patient_sv_name","patient_sv_merged","patient_sv_merged_coordinate","partner_sv_name","partner_sv_merged")) %>% 
    dplyr::mutate(sv_merged=patient_sv_merged,sv_merged_coordinate=patient_sv_merged_coordinate)
  
  #adjust the unfiltered_svs_df
  #might be confusing because the sv already seems linked? => no it is true, there was already a merged group that wasnt right.. otherwise not multitool
  #so should reset and reannotate
  unfiltered_svs_df_bk=unfiltered_svs_df
  unfiltered_svs_orphan = unfiltered_svs_df_bk %>% 
    filter(patient_sv_name %in% orphan_partner_svs_df$patient_sv_name) %>% 
    filter(partner_sv_name %in% orphan_source_svs_df$patient_sv_name)
  
  unfiltered_svs_orphan = unfiltered_svs_orphan %>% select(-patient_sv_merged,-sv_merged,-sv_merged_coordinate) %>% left_join(orphan_partner_svs_df)
  unfiltered_svs_orphan$flag_in_filtered_svs = T
  
  
  #also here adjust partner 
  unfiltered_svs_orphan_source = unfiltered_svs_df_bk %>% 
    filter(patient_sv_name %in% orphan_source_svs_df$patient_sv_name)
  unfiltered_svs_orphan_source = unfiltered_svs_orphan_source %>% select(-patient_sv_merged,-sv_merged,-sv_merged_coordinate) %>% left_join(orphan_source_svs_df %>% dplyr::mutate(sv_merged=patient_sv_merged,patient_sv_merged_coordinate=sv_merged_coordinate)
                                                                                                                                            %>% select(names(orphan_partner_svs_df)))
  unfiltered_svs_orphan_source$flag_in_filtered_svs = T
  
  #remove and add to dfs
  unfiltered_svs_df = unfiltered_svs_df_bk %>% 
    filter(!patient_sv_name %in% orphan_partner_svs_df$patient_sv_name) %>%
    filter(!patient_sv_name %in% orphan_source_svs_df$patient_sv_name)
  
  unfiltered_svs_df = rbind(unfiltered_svs_df,unfiltered_svs_orphan[,names(unfiltered_svs_df)],unfiltered_svs_orphan_source[,names(unfiltered_svs_df)])
  
  #make merged svs and add to df 
  orphan_partner_merged_svs = make_merged_svs(unfiltered_svs_orphan) %>% select(-partner_sv_merged,-partner_sv_merged_coordinate) %>% 
    left_join(orphan_partner_svs_df[,c("patient_sv_merged","partner_sv_merged")])
  
  orphan_partner_merged_svs = orphan_partner_merged_svs %>% left_join(orphan_merged_svs[,c("patient_sv_merged","sv_merged_coordinate")] %>% dplyr::rename(partner_sv_merged_coordinate=sv_merged_coordinate),by=c("partner_sv_merged"="patient_sv_merged")) %>% unique()
  
  orphan_partner_merged_svs = orphan_partner_merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  
  merged_svs = rbind(merged_svs,orphan_partner_merged_svs[,names(merged_svs)]) %>% unique()
  
  
  merged_svs %>% filter(svtype=="CTX" & !patient_sv_merged %in% merged_svs$partner_sv_merged) %>% nrow() ==0 
  merged_svs %>% filter(svtype=="CTX" & !partner_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() ==0 
  
  ## remove non autosomal
  ## TODO later look into ctx with non autosomes to see what we are filtering out
  non_autosomal_merged_svs = merged_svs %>% filter(!chrom %in% autosomes | (!is.na(partner_chrom) & partner_chrom!="" & !partner_chrom %in% autosomes))
  
  merged_svs = merged_svs %>% filter(chrom %in% autosomes & (is.na(partner_chrom) | partner_chrom=="" | partner_chrom %in% autosomes))
  
  #add partner coord for ctx matching
  svs_df = unfiltered_svs_df %>% filter(flag_in_filtered_svs) %>% filter(patient_sv_merged %in% merged_svs$patient_sv_merged)
  #check if merging goes wrong for orphans? => yes had to reannotate see above [fixed]
  svs_df = svs_df %>% left_join(merged_svs[,c("patient_sv_merged","sv_merged_coordinate","partner_sv_merged","partner_sv_merged_coordinate","cancer_type")])
  
  #check if these still are a thing:
  if(merged_svs %>% filter(svtype=="CTX") %>% filter(is.na(partner_sv_merged_coordinate)) %>% nrow() > 0) {
    print("Orphan svs still present?")
    print(merged_svs %>% filter(svtype=="CTX") %>% filter(is.na(partner_sv_merged_coordinate)) )
  }
  
  #check consistency merged svs and svs df 
  merged_svs %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% nrow() == 0
  unfiltered_svs_df %>% filter(!patient_sv_merged %in% svs_df$patient_sv_merged) %>% filter(patient_sv_merged %in% merged_svs$patient_sv_merged) %>% nrow() == 0
  
  if(as_global) {
    merged_svs <<- merged_svs
    unfiltered_svs_df <<- unfiltered_svs_df
    svs_df <<- svs_df
    non_autosomal_merged_svs <<- non_autosomal_merged_svs
  } else {
    print("Can only be used as global to override merged_svs, unfiltered_svs_df, svs_df")
  }
}
