parse_snv_geno_somatic = function(snv_vcf,run_tumor_id,run_normal_id) {
  snv_AF = as.data.frame(geno(snv_vcf)$AF)  %>% mutate(snv_id = rownames(.)) %>% 
    dplyr::rename("tumor_AF" := !!run_tumor_id) %>% dplyr::rename("normal_AF" := !!run_normal_id) 
  
  snv_AF = snv_AF %>% rowwise() %>% mutate(normal_AF = ifelse(length(normal_AF)>1,mean(as.double(unlist(normal_AF))),
                                                              as.double(normal_AF)),
                                           tumor_AF = ifelse(length(tumor_AF)>1,mean(as.double(unlist(tumor_AF))),
                                                             as.double(tumor_AF))) %>% ungroup()
  
  snv_DP = as.data.frame(geno(snv_vcf)$DP) %>% mutate(snv_id = rownames(.)) %>%
    dplyr::rename("tumor_DP" := !!run_tumor_id) %>% dplyr::rename("normal_DP" := !!run_normal_id)
  snv_AD = as.data.frame(geno(snv_vcf)$AD)  %>% mutate(snv_id = rownames(.)) %>%
    dplyr::rename("tumor_AD" := !!run_tumor_id) %>% dplyr::rename("normal_AD" := !!run_normal_id) %>%
    mutate(tumor_AD = as.character(tumor_AD),normal_AD = as.character(normal_AD))
  
  
  snv_geno_df = snv_AF %>%
    left_join(snv_AD,by="snv_id")  %>%
    left_join(snv_DP,by="snv_id") 
  
  return(snv_geno_df)
}
parse_snv_geno_germline = function(snv_vcf) {
  snv_AF_DP = as.data.frame(info(snv_vcf))  %>% mutate(snv_id = rownames(.)) %>% select(snv_id,AF,DP) %>%
    dplyr::rename(normal_AF =AF, normal_DP = DP) 
  
  #list length does not work on hpc and this is much faster
  snv_AF_DP[grepl(",",snv_AF_DP$normal_AF),] = snv_AF_DP[grepl(",",snv_AF_DP$normal_AF),] %>%
    rowwise() %>% mutate(normal_AF = mean(as.double(unlist(normal_AF)))) %>% ungroup() %>% as.data.frame()
  snv_AF_DP$normal_AF = as.double(snv_AF_DP$normal_AF)
  
  snv_geno_df = snv_AF_DP 
  
  return(snv_geno_df)
}

