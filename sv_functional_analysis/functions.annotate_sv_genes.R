## DEPRECATED??? 

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

## copied from match wgs functions but made safe for ensdb
annotate_metadata = function(all_gr,metadata_df) {
  all_gr=all_gr[order(names(all_gr))]
  metadata = mcols(all_gr,use.names = T)  
  metadata$bp_name = rownames(metadata)
  
  if("sourceId" %in% names(metadata_df) && nrow(dplyr::filter(metadata_df,sourceId %in% metadata$sourceId))>0) {
    metadata = metadata %>% merge(metadata_df %>% dplyr::select(sourceId,tumor_af,normal_af,somatic),by="sourceId",all.x=T) %>% unique()
  } 
  else if("bp_name" %in% names(metadata_df) && nrow(dplyr::filter(metadata_df,bp_name %in% metadata$bp_name))>0) {
    metadata = metadata %>% merge(metadata_df %>% dplyr::select(bp_name,tumor_af,normal_af,somatic),by="bp_name",all.x=T) %>% unique()
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
  
  somatic_vcf_geno$tumor_af = (as.integer(somatic_vcf_geno$tumor_VAR_PR) + as.integer(somatic_vcf_geno$tumor_VAR_SR)) / 
    ( as.integer(somatic_vcf_geno$tumor_REF_SR) + as.integer(somatic_vcf_geno$tumor_REF_PR) +
        as.integer(somatic_vcf_geno$tumor_VAR_SR) + as.integer(somatic_vcf_geno$tumor_VAR_PR) )
  
  somatic_vcf_geno$normal_af = (as.integer(somatic_vcf_geno$normal_VAR_PR) + as.integer(somatic_vcf_geno$normal_VAR_SR)) / 
    ( as.integer(somatic_vcf_geno$normal_REF_SR) + as.integer(somatic_vcf_geno$normal_REF_PR) +
        as.integer(somatic_vcf_geno$normal_VAR_SR) + as.integer(somatic_vcf_geno$normal_VAR_PR) )
  
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
  
  germline_vcf_geno$normal_af = (as.integer(germline_vcf_geno$normal_VAR_PR) + as.integer(germline_vcf_geno$normal_VAR_SR)) / 
    ( as.integer(germline_vcf_geno$normal_REF_SR) + as.integer(germline_vcf_geno$normal_REF_PR) +
        as.integer(germline_vcf_geno$normal_VAR_SR) + as.integer(germline_vcf_geno$normal_VAR_PR) )
  
  germline_vcf_geno$sourceId= gsub(".",":",germline_vcf_geno$sourceId,fixed=T)
  
  }
  
  ## Note: manta sometimes 0s for normal AF in diploid model, also somtimes normal > tumor AF in somatic 
  
  #breakpoints need to append and annotate with allele frequencies, tool, somatic or not 
  ##  filtering the sourceIds for "manta" to remove "svrecord" since they ofen overlap with normal variants, but without metadata
  
  somatic_gr = breakpointRanges(somatic_vcf)
  somatic_gr = somatic_gr[grepl("Manta",somatic_gr$sourceId),]
  somatic_gr = annotate_metadata(somatic_gr,somatic_vcf_geno)
  
  germline_gr = breakpointRanges(germline_vcf)
  germline_gr = germline_gr[grepl("Manta",germline_gr$sourceId),]
  germline_gr = annotate_metadata(germline_gr,germline_vcf_geno)
  
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
  gr_df = as.data.frame(all_gr) %>% mutate(bp_name = rownames(.)) %>%  dplyr::select(sourceId)
  
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
      as.data.frame(geno(vcf)$RV) %>% mutate(sourceId = rownames(.)) %>% 
        dplyr::rename("tumor_RV" := !!patient$tumor_id) %>% dplyr::rename("normal_RV" := !!patient$normal_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(vcf)$RR) %>% mutate(sourceId = rownames(.)) %>%
        dplyr::rename("tumor_RR" := !!patient$tumor_id) %>% dplyr::rename("normal_RR" := !!patient$normal_id)
      , by="sourceId"
    ) %>%
    mutate(tumor_af = (tumor_RV/(tumor_RR+tumor_RV))) %>%
    mutate(normal_af = (normal_RV/(normal_RR+normal_RV)))
  
  #if both are 0 then drop the NAs, they can be rescued by DVs or otherwise we just dont have data for that variant
  vcf_geno_precise = drop_na(vcf_geno_precise)
  
  # calculate based on DV for all to also rescue the NAs with RV
  vcf_geno_all = 
    vcf_geno %>% 
    left_join (
      as.data.frame(geno(vcf)$DV) %>% mutate(sourceId = rownames(.)) %>% 
        dplyr::rename("tumor_DV" := !!patient$tumor_id) %>% dplyr::rename("normal_DV" := !!patient$normal_id)
      , by="sourceId") %>%
    left_join(  
      as.data.frame(geno(vcf)$DR) %>% mutate(sourceId = rownames(.)) %>%
        dplyr::rename("tumor_DR" := !!patient$tumor_id) %>% dplyr::rename("normal_DR" := !!patient$normal_id)
      , by="sourceId"
    ) %>%
    mutate(tumor_af = (tumor_DV/(tumor_DR+tumor_DV))) %>%
    mutate(normal_af = (normal_DV/(normal_DR+normal_DV)))
  
  vcf_geno_all = drop_na(vcf_geno_all)
  
  #if AF in precise table, then exclude it here 
  vcf_af = rbind(vcf_geno_precise[,c("sourceId","tumor_af","normal_af")], 
                          vcf_geno_all[!vcf_geno_all$sourceId %in% vcf_geno_precise$sourceId,c("sourceId","tumor_af","normal_af")])
  
  
  ## DISABLED: If somatic file is provided, flag breakpoints as such by annotate metadata
  #input germline anno table, somatic file path, total intervals = GRanges() for delly 
  ##use the full DF for annotation here instead of the allele frequency one because of "drop NA"
  #gr_df = annotate_metadata_somatic(gr_df, argv$vcf_somatic, GRanges())
  
  ## Currently not used:
  gr_df$somatic = NA
  
  #merge allele frequencies, 
  gr_df = gr_df %>% left_join(vcf_af, by="sourceId")
  
  #add the annotation to breakpoints
  all_gr = annotate_metadata(all_gr,gr_df)
  
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
  gr_df = as.data.frame(all_gr) %>% mutate(bp_name = rownames(.))
  
  ## Allele frequency calculation
  #how to calculate AF depends on variant size 
  # svlen <1000 dont use the REFPAIR otherwise do 
  #uses the label of tumor/normal not the id.
  
  vcf_anno = gr_df %>% dplyr::select(bp_name,svLen) %>% left_join(
    as.data.frame(vcf_geno$VF) %>% mutate(bp_name = rownames(.)) %>% 
      dplyr::rename("tumor_VF" := !!patient$tumor_label) %>% dplyr::rename("normal_VF" := !!patient$normal_label)
    , by="bp_name") %>% left_join(
      as.data.frame(vcf_geno$REF) %>% mutate(bp_name = rownames(.)) %>% 
        dplyr::rename("tumor_REF" := !!patient$tumor_label) %>% dplyr::rename("normal_REF" := !!patient$normal_label)
      , by="bp_name") %>% left_join(
        as.data.frame(vcf_geno$REFPAIR) %>% mutate(bp_name = rownames(.)) %>% 
          dplyr::rename("tumor_REFPAIR" := !!patient$tumor_label) %>% dplyr::rename("normal_REFPAIR" := !!patient$normal_label)
        , by="bp_name") %>%
    mutate(tumor_af = ifelse(!is.na(svLen)&svLen<1000, (tumor_VF/(tumor_VF+tumor_REF)), (tumor_VF/(tumor_VF+tumor_REF+tumor_REFPAIR)))) %>%
    mutate(normal_af = ifelse(!is.na(svLen)&svLen<1000, (normal_VF/(normal_VF+normal_REF)), (normal_VF/(normal_VF+normal_REF+normal_REFPAIR))))
  
  ## DISABLED: annotate with somatic flag
  #vcf_anno = annotate_metadata_somatic(vcf_anno, argv$vcf_somatic, total_intervals)
  #rownames(vcf_anno) = vcf_anno$bp_name
  vcf_anno$somatic=NA
  
  #add the annotation to breakpoints
  all_gr = annotate_metadata(all_gr,vcf_anno)
  
  #annotate with tool
  mcols(all_gr)[["tool"]]="gridss"
  return(all_gr)
}


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

## deprecated

filter_gr_tumor_specific = function(gr) {
  gr=trim(gr)
  #prevent NAs
  gr_na_normal = gr[!is.na(gr$tumor_af)&is.na(gr$normal_af)]
  if(length(gr_na_normal)>0) {
    gr[!is.na(gr$tumor_af)&is.na(gr$normal_af)]$normal_af=0
  }
  gr$tumor_normal_af = gr$tumor_af-gr$normal_af
  gr$tumor_normal_af_ratio = gr$tumor_normal_af/gr$normal_af
  
  gr = gr[!is.na(gr$tumor_normal_af)&gr$tumor_normal_af>0.05&gr$tumor_normal_af_ratio>=1.5]
  
  return(gr)
}

## deprecated
anno_somatic_variant = function(summary_df) {
  summary_df = summary_df %>% mutate(tumor_normal_ratio = ((tumor_af-normal_af)/normal_af), 
                                     somatic_variant = (tumor_af-normal_af) > 0.05 & ((tumor_af-normal_af)/normal_af)>2)
  return(summary_df)
}


