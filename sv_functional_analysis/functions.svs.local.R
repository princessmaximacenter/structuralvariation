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

