#deprecated? replaced by other tryouts single category variable
#annotate sv af class 2021-08-19 
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


annotate_sv_af_class = function(df,tumor_normal_ratio_threshold=2) {
  # get rid of the flags and make mutually exclusive categories
  df = df %>% 
    dplyr::mutate(tumor_normal_diff = (tumor_af-normal_af)) %>% #not used
    dplyr::mutate(tumor_normal_ratio = (tumor_af/normal_af)) 
  
  remove_somatic_col=F
  if(!"somatic" %in% names(df)) {
    #add temp empty column for SV merged. Only SV unmerged has the 'somatic' column from Manta
    remove_somatic_col=T
    df$somatic=NA
  }
  
  #germline/somatic/low af flags are helper variables only 
  #todo: remove if dependancy on germline/somatic/low af flags is removed everywhere
  #use mantas somatic column if available 
  df = df %>% dplyr::mutate(germline_variant = (
                                                  (normal_af>=0.1 | (normal_af>0.05 & tumor_normal_ratio<tumor_normal_ratio_threshold))) |
                              (is.na(tumor_normal_ratio) & somatic==FALSE))
  
  #nongermline implies  normal_af<0.1 but not necessarily tumor_normal_ratio>tumor_normal_ratio_threshold 
  df = df %>% dplyr::mutate(somatic_variant = !germline_variant & 
                              ((!is.na(tumor_af) & tumor_af>0.05 & tumor_normal_ratio>=tumor_normal_ratio_threshold) | 
                                 (is.na(tumor_normal_ratio) & somatic==TRUE))) 
  
  df = df %>% dplyr::mutate(low_af = (!somatic_variant & !germline_variant & (tumor_af<0.1 & normal_af<0.1)) |
                              (tumor_af<0.1 & is.na(normal_af)) | 
                              (is.na(tumor_af)& normal_af<0.1) |
                              !(is.na(tumor_af)&is.na(normal_af)))
  
  df = df %>% mutate(variant_type = ifelse(somatic_variant==TRUE,"tumor_specific",
                                           ifelse(germline_variant==TRUE,"germline",
                                                  ifelse(low_af==TRUE,"low_af","ambiguous"))))
  df$variant_type = factor(df$variant_type)
  
  if(remove_somatic_col==T) {
    df = df %>% select(-somatic)
  }
  return(df)
}

get_variant_type_per_tool = function(svs_df) {
  #gives number of SVS per tool and the variant types associated with it
  
  svs_in_merged_tool = svs_df %>% 
    group_by(patient_sv_merged,tool,variant_type) %>%
    summarize(sv_cnt=n())
  return(svs_in_merged_tool)
}

get_variant_type_ambiguous_per_tool = function(svs_df) {
  ## are tools consistent or ambiguous, multiple variant types per sv merged?
  svs_in_merged_tool = get_variant_type_per_tool(svs_df)
  
  variant_type_ambiguous_per_tool = svs_in_merged_tool %>% ungroup() %>% 
    group_by(tool,patient_sv_merged) %>% 
    summarize(cnt = n(), variant_types=toString(unique(sort(variant_type)))) %>% filter(cnt>1) 
  return(variant_type_ambiguous_per_tool)
}

get_variant_type_majority_vote = function(svs_df) {
  svs_in_merged_tool = get_variant_type_per_tool(svs_df)
  
  #only count a tool once so delly with 5 snvs for germline is still 1 vote
  ## also see the ambiguity reporting function
  
  variant_type_majority_vote = svs_in_merged_tool %>% ungroup() %>% 
    group_by(patient_sv_merged) %>%
    count(patient_sv_merged, variant_type) %>% 
    summarise(variant_type = variant_type[n>sum(n/2)][1]) %>% as.data.frame()
  
  return(variant_type_majority_vote)
  
}


