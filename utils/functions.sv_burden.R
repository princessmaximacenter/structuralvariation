cnt_svs = function(svs_df, attr_name=NULL, group_cols=c("patient_id")) {
  cnt_table = svs_df %>% group_by(across(all_of(group_cols))) %>% summarise(svs_cnt = n(),.groups="keep") %>% ungroup()
  if(!is.null(attr_name)) {
    colnames(cnt_table) = c(group_cols,attr_name)
  }
  return(cnt_table)
}

get_sv_rows_per_tool = function(svs_df,group_cols=c("patient_id"),col_prefix="") {
  ##Note: just counts rows, so these are SV events for intra chr and breakpoints for CTX
  ## the select all cols except for tools makes it possible to count uq svs shared between tools
  if(nrow(svs_df)==0) { return(data.frame()) }
  
  patient_svs_cnt= cnt_svs(svs_df %>% dplyr::select(all_of(colnames(svs_df)[colnames(svs_df)!="tool"])) %>% unique(),attr_name=paste0(col_prefix,"all_cnt"),group_cols=group_cols) %>%
    left_join(cnt_svs(svs_df %>% filter(tool=="manta"),paste0(col_prefix,"manta_cnt"),group_cols=group_cols),by = group_cols) %>%
    left_join(cnt_svs(svs_df %>% filter(tool=="delly"),paste0(col_prefix,"delly_cnt"),group_cols=group_cols),by = group_cols) %>%
    left_join(cnt_svs(svs_df %>% filter(tool=="gridss"),paste0(col_prefix,"gridss_cnt"),group_cols=group_cols),by = group_cols) 
  
  patient_svs_cnt=as.data.frame(patient_svs_cnt)
  return(patient_svs_cnt)
}

## function get overall from sv type counts
sum_to_all_cnt_sv_type = function(sv_event_long) {
  if(nrow(sv_event_long)==0) { return(data.frame()) }
  
  sv_event_long=sv_event_long[sv_event_long$svtype!="all",]
  cols=names(sv_event_long)[!names(sv_event_long) %in% c("svtype","cnt")]
  sv_event_all=sv_event_long %>% group_by(across(all_of(cols))) %>% summarize(cnt=sum(cnt))
  sv_event_all$svtype="all"
  sv_event_all=as.data.frame(sv_event_all)
  return(sv_event_all)
}

convert_cnt_sv_row_breakpoint = function(svs_type_long) {
  #sv breakpoint is intra-chr *2 
  if(nrow(svs_type_long)==0) { return(data.frame()) }
  
  svs_type_long=svs_type_long[svs_type_long$svtype!="all",]
  sv_breakpoint_long = svs_type_long %>% mutate(cnt_type=str_replace(cnt_type,"_row","_breakpoint"),
                                                cnt= ifelse(svtype!="CTX",cnt*2,cnt))
  sv_breakpoint_long_all=sum_to_all_cnt_sv_type(sv_breakpoint_long)
  sv_breakpoint_long=rbind(as.data.frame(sv_breakpoint_long),as.data.frame(sv_breakpoint_long_all))
  return(sv_breakpoint_long)
}

convert_cnt_regional_sv_row_breakpoint = function(svs_type_orientation_long) {
  #sv overlaps with region using start/end coordinates 
  #sv orientation is taken into account => necessary because one of the SV bp can overlap a region only
  
  if(nrow(svs_type_orientation_long)==0) { return(data.frame()) }
  
  svs_type_orientation_long=svs_type_orientation_long[svs_type_orientation_long$svtype!="all",]
  
  #sum over the sv orientation column, adding start/end counts for each region
  cols=names(svs_type_orientation_long)[!names(svs_type_orientation_long) %in% c("sv_breakpoint_orientation","cnt")]
  svs_type_orientation_long=svs_type_orientation_long %>% group_by(across(all_of(cols))) %>% summarize(cnt=sum(cnt))
  
  #enforce start/stop also inflates CTX counts => should be /2 
  sv_breakpoint_long = svs_type_orientation_long %>% mutate(cnt_type=str_replace(cnt_type,"_row","_breakpoint"),
                                                            cnt= ifelse(svtype=="CTX",cnt/2,cnt))
  
  sv_breakpoint_long=as.data.frame(sv_breakpoint_long)
  sv_breakpoint_long_all=sum_to_all_cnt_sv_type(sv_breakpoint_long)
  sv_breakpoint_long=rbind(sv_breakpoint_long,as.data.frame(sv_breakpoint_long_all))
  
  return(sv_breakpoint_long)
}

convert_cnt_sv_row_event = function(svs_type_long) {
  #sv event is CTX/2
  if(nrow(svs_type_long)==0) { return(data.frame()) }
  svs_type_long=svs_type_long[svs_type_long$svtype!="all",]
  sv_event_long = svs_type_long %>% mutate(cnt_type=str_replace(cnt_type,"_row","_event"),
                                           cnt= ifelse(svtype=="CTX",cnt/2,cnt))
  sv_event_long=as.data.frame(sv_event_long)
  sv_event_long_all=sum_to_all_cnt_sv_type(sv_event_long)
  sv_event_long=rbind(sv_event_long,as.data.frame(sv_event_long_all))
  
  return(sv_event_long)
}

get_sv_burden_long = function(svs_df,group_cols=c("patient_id"),cnt_type="sv_row") {
  #long datafame 
  
  if(nrow(svs_df)==0) {return (data.frame())}
  svs_long=get_sv_rows_per_tool(svs_df,group_cols = group_cols)
  svs_long[is.na(svs_long)]=0 #convert NAs to 0 counts
  if(!"svtype" %in% group_cols) svs_long$svtype="all"
  svs_long$cnt_type=cnt_type
  svs_long = svs_long %>% pivot_longer(contains("_cnt"), names_to="tool",values_to="cnt")
  
  svs_long_all=sum_to_all_cnt_sv_type(svs_long)
  sv_rows_long=rbind(as.data.frame(svs_long_all),as.data.frame(svs_long))
  return(sv_rows_long)
}

get_merged_sv_binned_event_long = function(multitool_svs_df) {
  
  
  ## bin by sv size
  merged_sv_len_df = multitool_svs_df %>% group_by(patient_sv_merged,svtype) %>% summarize(svLen=mean(svLen,na.rm=T))
  multitool_svs_merged_svlen_df = multitool_svs_merged_df %>% left_join(merged_sv_len_df)
  multitool_svs_merged_svlen_df = make_svlen_bin(multitool_svs_merged_svlen_df)
  
  binned_group_cols=c(group_cols,"svlen_bin")
  patient_merged_sv_binned_svlen_cnt_long = get_sv_burden_long(multitool_svs_merged_svlen_df,group_cols=binned_group_cols,cnt_type="merged_row")
  merged_sv_binned_svlen_event_long = patient_merged_sv_binned_svlen_cnt_long %>% convert_cnt_sv_row_event()
  
  #check merged_sv_binned_svlen_event_long %>% group_by(patient_id,cnt_type,tool,svtype) %>% summarize(cnt = sum(cnt)) %>% arrange(svtype,tool) == merged_sv_event_long%>% arrange(svtype,tool)
  #for union merged_sv_binned_svlen_event_long %>% group_by(patient_id,cnt_type,tool,svtype,variant_type) %>% summarize(cnt = sum(cnt)) %>% arrange(svtype,tool)
  
  ## bin by multi sv tool support
  merged_sv_tool_df = multitool_svs_df %>% group_by(patient_sv_merged) %>% summarize(tool_bin=toString(sort(unique(tool))))
  multitool_svs_merged_tool_df = multitool_svs_merged_df %>% left_join(merged_sv_tool_df)
  binned_group_cols=c(group_cols,"tool_bin")
  patient_merged_sv_binned_tool_cnt_long = get_sv_burden_long(multitool_svs_merged_tool_df,
                                                              group_cols=binned_group_cols,cnt_type="merged_row")
  
  merged_sv_binned_tool_event_long = patient_merged_sv_binned_tool_cnt_long %>% convert_cnt_sv_row_event() %>% filter(tool=="all_cnt")
  #checked with the upset plot /venn and is correct 
  
  ## bin by VAF 
  merged_sv_vaf_df = multitool_svs_df %>% group_by(patient_sv_merged,svtype) %>% summarize(tumor_af=mean(tumor_af,na.rm=T))
  
  multitool_svs_merged_vaf_df = multitool_svs_merged_df %>% left_join(merged_sv_vaf_df)
  multitool_svs_merged_vaf_df = make_tumor_af_bin(multitool_svs_merged_vaf_df)
  
  binned_group_cols=c(group_cols,"tumor_af_bin")
  patient_merged_sv_binned_vaf_cnt_long = get_sv_burden_long(multitool_svs_merged_vaf_df,group_cols=binned_group_cols,cnt_type="merged_row")
  merged_sv_binned_vaf_event_long = patient_merged_sv_binned_vaf_cnt_long %>% convert_cnt_sv_row_event()
  
  #check 
  #merged_sv_binned_vaf_event_long %>% group_by(patient_id,cnt_type,tool,svtype) %>% summarize(cnt = sum(cnt)) %>% arrange(svtype,tool) == merged_sv_event_long%>% arrange(svtype,tool)
  #merged_sv_binned_vaf_event_long %>% filter(tumor_af_bin=="0.1<taf<0.2" &svtype!="all")
  #merged_sv_vaf_df %>% dplyr::select(patient_sv_merged,svtype,tumor_af) %>% unique() %>% filter(tumor_af>=0.1&tumor_af<0.2) %>% group_by(svtype) %>% summarize(n())
  
  #check for variant type binning 
  #merged_sv_vaf_df %>% dplyr::select(patient_sv_merged,svtype,tumor_af) %>% unique() %>% left_join(svs_df[,c("patient_sv_merged","variant_type")]) %>% unique() %>% filter(tumor_af>=0.1&tumor_af<0.2) %>% group_by(svtype,variant_type) %>% summarize(n())
  #pivot longer 
  #bin_key, bin_value 
  merged_sv_binned_svlen_event_long=merged_sv_binned_svlen_event_long %>% pivot_longer(svlen_bin,names_to="bin_key",values_to="bin_value")
  merged_sv_binned_tool_event_long=merged_sv_binned_tool_event_long %>% pivot_longer(tool_bin,names_to="bin_key",values_to="bin_value")
  merged_sv_binned_vaf_event_long=merged_sv_binned_vaf_event_long %>% pivot_longer(tumor_af_bin,names_to="bin_key",values_to="bin_value")
  
  merged_sv_binned_event_long = rbind(merged_sv_binned_svlen_event_long,merged_sv_binned_tool_event_long,merged_sv_binned_vaf_event_long)
  return(merged_sv_binned_event_long)
}


get_sv_binned_event_long = function(svs_df) {
  #for all svs
  binned_svs_df = svs_df %>% make_svlen_bin() %>% make_tumor_af_bin()
  
  patient_sv_binned_svlen_cnt_long = get_sv_burden_long(binned_svs_df,group_cols=c(group_cols,"svlen_bin"),cnt_type="sv_row")
  sv_binned_svlen_event_long = patient_sv_binned_svlen_cnt_long %>% convert_cnt_sv_row_event()
  
  patient_sv_binned_vaf_cnt_long = get_sv_burden_long(binned_svs_df,group_cols=c(group_cols,"tumor_af_bin"),cnt_type="sv_row")
  sv_binned_vaf_event_long = patient_sv_binned_vaf_cnt_long %>% convert_cnt_sv_row_event()
  
  #sv_binned_vaf_event_long %>% filter(tumor_af_bin=="0.1<taf<0.2" &svtype!="all")
  #svs_df %>% dplyr::select(patient_sv_name,svtype,tumor_af)  %>% filter(tumor_af>=0.1&tumor_af<0.2) %>% group_by(svtype) %>% summarize(n())
  
  #note that both gridss bp can have different AF for CTX
  #svs_df %>% dplyr::select(patient_sv_name,svtype,tumor_af)  %>% filter(tumor_af<0.1 & svtype=="CTX")
  
  #pivot longer 
  #bin_key, bin_value 
  sv_binned_svlen_event_long=sv_binned_svlen_event_long %>% pivot_longer(svlen_bin,names_to="bin_key",values_to="bin_value")
  sv_binned_vaf_event_long=sv_binned_vaf_event_long %>% pivot_longer(tumor_af_bin,names_to="bin_key",values_to="bin_value")
  
  sv_binned_event_long = rbind(sv_binned_svlen_event_long,sv_binned_vaf_event_long)
  return(sv_binned_event_long)
}
