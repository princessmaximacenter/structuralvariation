#Note at december 2021 not sure if still used
load_cna_seg = function(cna_path) {
  cna_file = Sys.glob(cna_path)
  if(length(cna_file)!=1){ 
    print("WARNING: input unclear")
    print(cna_file)
    return(GRanges)
  }
  
  segments = read.table(cna_file,sep="\t",comment.char="@",header=T,skip = 1)
  segments=segments %>% dplyr::rename(seqnames = CONTIG, start=START, end=END,copy_ratio_log2_mean=MEAN_LOG2_COPY_RATIO,call=CALL)
  cna = GRanges(segments) 
  names(cna)=paste0("cna_",1:length(cna))
  mcols(cna)$coordinate = paste0(seqnames(cna),":",start(cna),"-",end(cna))
  return(cna)
}

#annotate cna chr and gene
#december 2021
read_modeled_seg = function(segments_file){
  modeled_seg = read.table(segments_file,comment.char = "@",header=T,sep="\t")
  modeled_seg = modeled_seg %>% 
    dplyr::rename(seqnames=CONTIG,start=START,end=END,
                  maf_50=MINOR_ALLELE_FRACTION_POSTERIOR_50,cr_l2fc_50=LOG2_COPY_RATIO_POSTERIOR_50,
                  maf_10=MINOR_ALLELE_FRACTION_POSTERIOR_10,cr_l2fc_10=LOG2_COPY_RATIO_POSTERIOR_10,
                  maf_90=MINOR_ALLELE_FRACTION_POSTERIOR_90,cr_l2fc_90=LOG2_COPY_RATIO_POSTERIOR_90)
  
  
  modeled_seg$cna_id=paste0("cna_",1:nrow(modeled_seg))
  modeled_seg$width=modeled_seg$end-modeled_seg$start
  return(modeled_seg)
}

## default
get_overlaps_modeled_seg_summary = function(overlaps_df,group_cols) {
  summary = overlaps_df %>% group_by(across(all_of(group_cols))) %>%
    summarize(cr_l2fc_50_max = max(cr_l2fc_50,na.rm = T),
              maf_50 = sum(maf_50*overlap_set2_set1,na.rm = T), 
              cr_l2fc_50 = sum(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              maf_10 = sum(maf_10*overlap_set2_set1,na.rm = T), 
              cr_l2fc_10 = sum(cr_l2fc_10*overlap_set2_set1,na.rm = T),
              maf_90 = sum(maf_90*overlap_set2_set1,na.rm = T), 
              cr_l2fc_90 = sum(cr_l2fc_90*overlap_set2_set1,na.rm = T),
              seg_cnt=n(), loss_cnt=sum(call=="loss"), gain_cnt=sum(call=="gain"), loh_cnt=sum(call=="loh"),
              seg_covered = sum(overlap_set2_set1),
              .groups="drop") 
  summary = summary %>% as.data.frame()
  
  ## Note that : if MAF is NA everywhere then becomes 0 and then later on causes LOH call. -> reset afterwards to not mess up the function
  check_missing_maf_data = overlaps_df %>% 
    filter(!is.na(maf_50)) %>% 
    group_by(across(all_of(group_cols))) %>% summarize(seg_covered_maf_data = sum(overlap_set2_set1))
  
  check_missing_cr_data = overlaps_df %>% 
    filter(!is.na(cr_l2fc_50)) %>% 
    group_by(across(all_of(group_cols))) %>% summarize(seg_covered_cr_data = sum(overlap_set2_set1))
  
  summary = summary %>% left_join(check_missing_maf_data) %>% left_join(check_missing_cr_data)
  summary[is.na(summary$seg_covered_maf_data),c("seg_covered_maf_data")]=0
  summary[is.na(summary$seg_covered_cr_data),c("seg_covered_cr_data")]=0
  
  summary[summary$seg_covered_maf_data==0,c("maf_50","maf_90","maf_10")] = NA
  
  return(summary)

}

## short
get_overlaps_modeled_seg_summary_short = function(overlaps_df,group_cols) {
  summary = overlaps_df %>% group_by(across(all_of(group_cols))) %>%
    summarize(cr_l2fc_50_max = max(cr_l2fc_50,na.rm = T),
              maf_50 = sum(maf_50*overlap_set2_set1,na.rm = T), 
              cr_l2fc_50 = sum(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              # cr_sd = sd(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              # cr_var = var(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              # cr_iqr = iqr(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              
              seg_cnt=n(), 
              seg_covered = sum(overlap_set2_set1),
              .groups="drop") 
  summary = summary %>% as.data.frame()
  
  ## Note that : if MAF is NA everywhere then becomes 0 and then later on causes LOH call. -> reset afterwards to not mess up the function
  check_missing_maf_data = overlaps_df %>% 
    filter(!is.na(maf_50)) %>% 
    group_by(across(all_of(group_cols))) %>% summarize(seg_covered_maf_data = sum(overlap_set2_set1))
  
  check_missing_cr_data = overlaps_df %>% 
    filter(!is.na(cr_l2fc_50)) %>% 
    group_by(across(all_of(group_cols))) %>% summarize(seg_covered_cr_data = sum(overlap_set2_set1))
  
  summary = summary %>% left_join(check_missing_maf_data) %>% left_join(check_missing_cr_data)
  summary[is.na(summary$seg_covered_maf_data),c("seg_covered_maf_data")]=0
  summary[is.na(summary$seg_covered_cr_data),c("seg_covered_cr_data")]=0
  
  summary[summary$seg_covered_maf_data==0,c("maf_50")] = NA
  
  return(summary)
  
}

get_relative_cr_l2fc = function(cohort_chr_arm_cna,cohort, cohort_col=c("tumor_id","cr_mean_weighted","ploidy_status")){
  if(!(all(cohort_col %in% names(cohort)))) {
    print(paste0("Missing columns cohort df, need ",toString(cohort_col)))
    return(cohort_chr_arm_cna)
  }
  cohort_chr_arm_cna_relative = cohort_chr_arm_cna %>% 
    left_join(cohort[,cohort_col]) %>%
    dplyr::rename(call_absolute=call)%>% #safeguard
    mutate(relative_cr_l2fc_50=ifelse(ploidy_status=="poly",cr_l2fc_50-cr_mean_weighted,cr_l2fc_50)) %>% 
    call_cna(cna_colname="relative_cr_l2fc_50",maf_colname = "maf_50") %>%
    dplyr::rename(relative_call=call, call=call_absolute) 
  
  return(cohort_chr_arm_cna_relative)
}

rename_relative_cr_to_normal = function(df) {
  cols = c("call","cr_l2fc_50","relative_call","relative_cr_l2fc_50")
  if(length(cols[!cols %in% colnames(df)]) == 0) {
    df = df %>% select(-call,-cr_l2fc_50) %>% dplyr::rename(call=relative_call,cr_l2fc_50=relative_cr_l2fc_50)
  } else {
    print("cols not present, did not rename")
  }
  return(df)
}

get_contig_lengths  =  function(segments_file,chromosomes=c(paste("chr",1:22,sep=""))){
  contig_lengths = data.frame()
  for (line in readLines(segments_file,)){
    if(startsWith(line,"@SQ")) {
      columns = strsplit(line,"\t")[[1]]
      chromosome = columns[grepl("SN:",columns)]
      length = columns[grepl("LN:",columns)]
      
      contig=c()
      contig$seqnames = str_replace(chromosome,"SN:","")
      contig$length = str_replace(length,"LN:","")
      
      #by default skip chrX,chrY, alt, chrUn and decoy 
      if(contig$seqnames %in% chromosomes) { 
        contig_lengths = rbind(contig_lengths,contig)
      } else { next() }
    } else if (startsWith(line,"chr")){
      #end of header
      break
    }
  }
  contig_lengths$length=as.numeric(contig_lengths$length)
  
  return(contig_lengths)
}

get_fga_per_chrom = function(modeled_seg,contig_lengths,cr_threshold=0.2,flag_use_provided_length=F) {
  if(flag_use_provided_length==T) {
    print("Using lengths of input segments instead of contigs")
    contig_lengths=modeled_seg %>% group_by(seqnames) %>% summarize(length=sum(width, na.rm=T)) %>% as.data.frame()
  }
  fga_per_chrom = data.frame()
  for(chrom in contig_lengths$seqnames) {
    cna_segments = modeled_seg[modeled_seg$seqnames==chrom & abs(modeled_seg$cr_l2fc_50) > cr_threshold,]
    
    fga_entry = c()
    fga_entry$seqnames = chrom
    fga_entry$total_bases = sum(modeled_seg[modeled_seg$seqnames==chrom,]$width,na.rm = T)
    fga_entry$cna_bases = sum(cna_segments$width)
    fga_entry$seg_cnt = nrow(modeled_seg[modeled_seg$seqnames==chrom,])
    fga_entry$seg_cna_cnt = nrow(cna_segments)
    fga_entry$fga = fga_entry$cna_bases / contig_lengths[contig_lengths$seqnames==chrom,c("length")]
    fga_entry$seg_covered = fga_entry$total_bases / contig_lengths[contig_lengths$seqnames==chrom,c("length")]
    fga_entry$max_cna = ifelse(fga_entry$cna_bases>0,max(cna_segments$cr_l2fc_50),0)
    fga_entry$min_cna = ifelse(fga_entry$cna_bases>0,min(cna_segments$cr_l2fc_50),0)
    
    fga_per_chrom = rbind(fga_per_chrom,fga_entry)
    
  }
  
  ## Overall FGA
  fga_entry = c()
  fga_entry$seqnames = "total"
  fga_entry$total_bases = sum(modeled_seg$width,na.rm = T)
  fga_entry$cna_bases = sum(fga_per_chrom$cna_bases)
  fga_entry$seg_cnt = sum(fga_per_chrom$seg_cnt)
  fga_entry$seg_cna_cnt = sum(fga_per_chrom$seg_cna_cnt)
  fga_entry$fga = fga_entry$cna_bases / sum(contig_lengths$length)
  fga_entry$seg_covered = fga_entry$total_bases / sum(contig_lengths$length)
  fga_entry$max_cna = max(fga_per_chrom$max_cna)
  fga_entry$min_cna = min(fga_per_chrom$min_cna)
  
  fga_per_chrom = rbind(fga_per_chrom,fga_entry)
  
  return(fga_per_chrom)
}
get_fga_per_chrom_loh = function(modeled_seg,contig_lengths,maf_threshold=0.4,flag_use_provided_length=F) {
  if(flag_use_provided_length==T) {
    print("Using lengths of input segments instead of contigs")
    contig_lengths=modeled_seg %>% group_by(seqnames) %>% summarize(length=sum(width, na.rm=T)) %>% as.data.frame()
  }
  fga_per_chrom = data.frame()
  modeled_seg=modeled_seg[!is.na(modeled_seg$maf_50),]
  for(chrom in contig_lengths$seqnames) {
    cna_segments = modeled_seg[modeled_seg$seqnames==chrom & modeled_seg$maf_50 < maf_threshold,]
    
    fga_entry = c()
    fga_entry$seqnames = chrom
    fga_entry$total_bases = sum(modeled_seg[modeled_seg$seqnames==chrom,]$width,na.rm = T)
    fga_entry$loh_bases = sum(cna_segments$width,na.rm = T)
    fga_entry$fga_loh = fga_entry$loh_bases / contig_lengths[contig_lengths$seqnames==chrom,c("length")]
    fga_entry$seg_covered = fga_entry$total_bases / contig_lengths[contig_lengths$seqnames==chrom,c("length")]
    
    fga_per_chrom = rbind(fga_per_chrom,fga_entry)
    
  }
  
  ## Overall FGA
  fga_entry = c()
  fga_entry$seqnames = "total"
  fga_entry$total_bases = sum(modeled_seg$width,na.rm = T)
  fga_entry$loh_bases = sum(fga_per_chrom$loh_bases,na.rm = T)
  fga_entry$fga_loh = fga_entry$loh_bases / sum(contig_lengths$length)
  fga_entry$seg_covered = fga_entry$total_bases / sum(contig_lengths$length)
  fga_per_chrom = rbind(fga_per_chrom,fga_entry)
  
  return(fga_per_chrom)
}

call_cna = function(cna, cna_colname="cr_l2fc_50", cna_cutoff=0.2, maf_colname="maf_50", maf_cutoff=0.4) {
  convert_gr = F
  if(isS4(cna)) { 
    cna_names = names(cna)
    cna = as.data.frame(cna) 
    convert_gr=T
  }
  cna$call=NA
  if(!is.null(cna_colname) & cna_colname %in% names(cna)) {
    cna[!is.na(cna[,cna_colname]) & cna[,cna_colname]>cna_cutoff,c("call")] = "gain"
    cna[!is.na(cna[,cna_colname]) & cna[,cna_colname]<(-cna_cutoff),c("call")] = "loss"
  }
  if(!is.null(maf_colname) & maf_colname %in% names(cna)) {
    cna[is.na(cna$call) & !is.na(cna[,maf_colname]) &  cna[,maf_colname]<maf_cutoff,c("call")] = "loh"
  }
  cna[is.na(cna$call),c("call")] = "0"
  
  if(convert_gr) {
    cna = GRanges(cna)
    names(cna)=cna_names
  }
  return(cna)
}

call_cn_seg_same_state = function(overlaps_df,summary_df,group_cols="chr_arm",cr_deviation_threshold=0.33 ) {
  #adapted from chr arms
  #stability by percentage of bases 
  ## cr threshold 1.33 and 0.67, for MAF only max because often jumping around between 0 and 0.5 
  ## how to do this for neutral chromosomes with bit of noise? => the min/max 0.1
  ## mean_cr_l2fc_50 is misleading if low seg covered but acen/stalk/gvar excluded anyway
  
  thresh_cr = summary_df %>% 
    dplyr::rename(mean_cr_l2fc_50=cr_l2fc_50) %>% 
    dplyr::select(group_cols,mean_cr_l2fc_50) %>%
    mutate(min_cr = ifelse(abs(mean_cr_l2fc_50)<0.1, -0.1, (1-cr_deviation_threshold)*mean_cr_l2fc_50),
           max_cr = ifelse(abs(mean_cr_l2fc_50)<0.1, 0.1, (1+cr_deviation_threshold)*mean_cr_l2fc_50))
  
  #min max should be opposite if negative numbers 
  chr_arm_cr_stability = overlaps_df %>% left_join(thresh_cr,by=group_cols) %>%
    group_by(across(all_of(c(group_cols,"mean_cr_l2fc_50")))) %>%
    filter( 
      (min_cr<max_cr & cr_l2fc_50 > min_cr & cr_l2fc_50 < max_cr) | 
        (max_cr<min_cr & cr_l2fc_50 < min_cr & cr_l2fc_50 > max_cr)
    ) %>%
    summarize(seg_covered_cr_stable = sum(overlap_set2_set1))
  
  summary_df = summary_df  %>% 
    left_join(chr_arm_cr_stability[,c(group_cols,"seg_covered_cr_stable")],by=group_cols)
  
  summary_df[is.na(summary_df$seg_covered_cr_stable),c("seg_covered_cr_stable")]=0
  
  if("maf_50" %in% names(summary_df)) {
    thresh_maf = summary_df %>% 
      dplyr::rename(mean_maf_50=maf_50) %>% 
      dplyr::select(group_cols,mean_maf_50) %>%
      mutate(# min_maf = (1-(cr_deviation_threshold/2))*mean_maf_50,
        max_maf = max(0.1,(1+(cr_deviation_threshold/2))*mean_maf_50))
    
    chr_arm_maf_stability = overlaps_df %>% left_join(thresh_maf,by=group_cols) %>%
      group_by(across(all_of(c(group_cols,"mean_maf_50")))) %>%
      filter(  maf_50 < max_maf) %>%
      summarize(seg_covered_maf_stable = sum(overlap_set2_set1))
    
    summary_df = summary_df  %>% 
      left_join(chr_arm_maf_stability[,c(group_cols,"seg_covered_maf_stable")],by=group_cols)
    summary_df[is.na(summary_df$seg_covered_maf_stable),c("seg_covered_maf_stable")]=0
    
  }
  
  
  return(summary_df)
}


call_cn_range_fraction = function(overlaps_df,range_length_col="width",group_cols=c("patient_label"),flag_loh_as_neutral=F) {
  overlaps_df=overlaps_df %>% unique() #needed!
  
  if(flag_loh_as_neutral) {
    overlaps_df[overlaps_df$call=="loh",c("call")]="0"
  }
  #weighted average sum(weight*value)/sum(weights)
  #works if I then use sum(cr_mean*frac_covered) afterwards over call states to get chr arm level mean
  summary = overlaps_df %>% group_by(across(all_of(c("call",group_cols)))) %>%
    summarize(seg_cnt=n(), 
              frac_covered = sum(overlap_set2_set1),
              seg_covered := sum(overlap_set2_set1*!!sym(range_length_col)),
              cr_l2fc_50 = sum(overlap_set2_set1*cr_l2fc_50)/sum(overlap_set2_set1))
  return(summary)
}

get_cn_status_windows = function(windows_df,cnas,modeled_seg,ranges_id_col="seqnames",windows_df_cols=NULL) {
  if(is.null(windows_df_cols)) {
    windows_df_cols=c("patient_sv_merged","window_orientation","window_id","window_coordinate","window_width",ranges_id_col) 
  }
  
  if(length(windows_df_cols[!windows_df_cols %in% names(windows_df)])>0){
    print("Missing cols windows df")
    print(toString(windows_df_cols[!windows_df_cols %in% names(windows_df)]))
    return(data.frame())
  }
  windows_df = windows_df[,c("seqnames","start","end",windows_df_cols)] %>% unique()
  windows=GRanges(windows_df)
  names(windows)=windows_df$window_id
  
  
  overlaps_cna_window = get_reciprocal_overlap_pairs(cnas,windows,reciprocal_overlap = 0,svtype_matching = F)
  overlaps_cna_window = overlaps_cna_window %>% 
    dplyr::rename(window_id=set2,cna_id=set1) %>% 
    left_join(modeled_seg[,c(cnas_df_cols)])  %>% 
    left_join(windows_df[,c(windows_df_cols)]) 
  
  cna_window_group_cols=windows_df_cols
  
  #todo function get_overlaps_modeled_seg_summary_short 
  cna_windows = overlaps_cna_window %>% group_by(across(all_of(cna_window_group_cols))) %>%
    summarize(cr_l2fc_50_max = max(cr_l2fc_50,na.rm = T),
              maf_50 = sum(maf_50*overlap_set2_set1,na.rm = T), 
              cr_l2fc_50 = sum(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              cr_sd = sd(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              seg_cnt=n(), 
              loss_cnt=sum(call=="loss"), gain_cnt=sum(call=="gain"), loh_cnt=sum(call=="loh"),
              seg_covered = sum(overlap_set2_set1),
              .groups="drop") 
  cna_windows = cna_windows  %>% call_cna()
  
  cna_windows = call_cn_seg_same_state(overlaps_cna_window,cna_windows,group_cols=c("window_id")) %>%
    call_cn_stability() 
  #deprecated use below instead %>% mutate(window_in_cr_state = (seg_covered_cr_stable*window_width))
  
  #get fraction of window gain/loss/0 ignore LOH
  cna_windows_cn_fractions = call_cn_range_fraction(overlaps_cna_window,range_length_col = "window_width",group_cols=c("window_id"),flag_loh_as_neutral=T)
  
  #need to divide by fraction to know whether all gain segments are within certain value relative to all gains instead of entire window.
  cna_windows_cn_fractions = call_cn_seg_same_state(overlaps_cna_window,cna_windows_cn_fractions,group_cols=c("window_id","call")) %>%
    dplyr::mutate(call_seg_covered_cr_stable=seg_covered_cr_stable/frac_covered)
  
  cna_windows_cn_fractions_wide = cna_windows_cn_fractions %>% 
    dplyr::select(window_id,call,frac_covered,cr_l2fc_50,call_seg_covered_cr_stable) %>% 
    pivot_wider(names_from="call",values_from=c("frac_covered","cr_l2fc_50","call_seg_covered_cr_stable")) 
  
  cna_windows = cna_windows %>% left_join(cna_windows_cn_fractions_wide)
  return(cna_windows)
}

call_cn_stability = function(df) {
  df[is.na(df$seg_covered_cr_stable),c("seg_covered_cr_stable")]=0
  df[is.na(df$seg_covered_maf_stable),c("seg_covered_maf_stable")]=0
  df = df %>% mutate(cr_stable = seg_covered_cr_stable>0.7,
                     maf_stable = seg_covered_maf_stable>0.6)
  return(df)
}

call_cn_alteration = function(df) {
  df = df %>% mutate(
    alteration = ifelse(cr_l2fc_50>0.2 & cr_stable,"gain_stable",
                        ifelse(cr_l2fc_50<(-0.2) & cr_stable,"loss_stable",
                               ifelse(cr_l2fc_50>0.2 & !cr_stable,"gain_unstable",
                                      ifelse(cr_l2fc_50<(-0.2) & !cr_stable,"loss_unstable",
                                             ifelse(abs(cr_l2fc_50)<0.2 & maf_50<0.4 & cr_stable, "cn_loh",
                                                    ifelse(abs(cr_l2fc_50)<0.2 & cr_stable, "0","unknown")))))))
  
  return(df)
}

call_cn_copies = function(df) {
  if(!"alteration" %in% names(df)) {
    df = df %>% call_cn_stability() %>% call_cn_alteration()
  }
  #determine min/max per patient with cap at expect when 1 copy gained/lost at 100% tumor cell percentage => used for scaling
  #log2(1.5) ~0.58
  #log2(0.5) = -1
  
  #only take into account stable chrom 
  summary_relative_cn = df %>% filter(alteration %in% c("gain_stable","loss_stable")) %>% 
    group_by(patient_id) %>%
    summarize(max_cr=min(0.6,max(cr_l2fc_50,na.rm = T)),
              min_cr=max(-1,min(cr_l2fc_50,na.rm = T)),
              min_maf=min(maf_50,na.rm = T))
  summary_relative_cn = add_missing_patients_to_overview(summary_relative_cn,df)
  summary_relative_cn[is.na(summary_relative_cn$max_cr),c("max_cr")]=0.6
  summary_relative_cn[is.na(summary_relative_cn$min_cr),c("min_cr")]=-1
  
  df = df %>% left_join(summary_relative_cn)
  
  df = df %>% mutate(
    nr_of_copies = round(2^cr_l2fc_50,1)*2,
    relative_cr = ifelse(grepl("gain",alteration),(cr_l2fc_50/max_cr),ifelse(grepl("loss",alteration),(cr_l2fc_50/min_cr),NA)),
    heterozygous_expected = ifelse(grepl("gain",alteration),log2(1.5),ifelse(grepl("loss",alteration),-1,NA)),
    relative_copies = round(2^(relative_cr*heterozygous_expected),1)*2,
    tumor_af = ifelse(grepl("_stable",alteration),round(abs((cr_l2fc_50/heterozygous_expected)/2),2),NA))
  
  #  tumor_af_adj = ifelse(grepl("_stable",alteration),round(abs(relative_cr)/2,2),NA))
  
  remove_cols = c(names(summary_relative_cn),"heterozygous_expected")
  remove_cols = remove_cols[!remove_cols %in% c("patient_id")]
  df = df %>% dplyr::select(-all_of(remove_cols))
  return(df)
  
}

get_sv_cn_average = function(svs,cnas,svs_df,cnas_df,svs_df_cols=NULL,cnas_df_cols=NULL){
  if(is.null(svs_df_cols)) {
    svs_df_overlap_cols = svs_df_overlap_cols[svs_df_overlap_cols %in% names(svs_df)]
    svs_df_anno_cols = svs_df_anno_cols[svs_df_anno_cols %in% names(svs_df)]
    svs_df_cols=c(svs_df_overlap_cols,svs_df_anno_cols)
  } 
  if(is.null(cnas_df_cols)) {
    cnas_df_cols = modeled_seg_cols[modeled_seg_cols %in% names(cnas_df)]
  }
  
  #because get_overlaps_modeled_seg_summary uses overlap_set2_set1 by default, have cnas as first 
  sv_cna_overlaps= get_reciprocal_overlap_pairs(cnas,svs,reciprocal_overlap = 0,svtype_matching = FALSE)
  
  if(nrow(sv_cna_overlaps)==0) { return(data.frame()) }
  join_cna_df = cnas_df[,c(cnas_df_cols,"start","end")] %>% dplyr::rename(cna_start=start,cna_end=end) %>% unique()
  join_svs_df = svs_df[,c(svs_df_cols,"start","end")] %>% dplyr::rename(sv_start=start,sv_end=end) %>% unique()
  
  sv_cna_overlaps = sv_cna_overlaps %>% dplyr::rename(patient_sv_merged=set2, cna_id=set1) %>% 
    left_join(join_cna_df)  %>%
    left_join(join_svs_df) 
  sv_cna = get_overlaps_modeled_seg_summary_short(sv_cna_overlaps,group_cols=c("patient_sv_merged"))
  sv_cna = sv_cna %>% call_cna()
  
  sv_cna = call_cn_seg_same_state(sv_cna_overlaps,sv_cna,group_cols=c("patient_sv_merged")) %>%
    call_cn_stability()
  
  #get fraction of window gain/loss/0 ignore LOH
  sv_cna_fractions = call_cn_range_fraction(sv_cna_overlaps,range_length_col = "svlen_mean",group_cols=c("patient_sv_merged"),flag_loh_as_neutral=T)
  
  #need to divide by fraction to know whether all gain segments are within certain value relative to all gains instead of entire window.
  sv_cna_fractions = call_cn_seg_same_state(sv_cna_overlaps,sv_cna_fractions,group_cols=c("patient_sv_merged","call")) %>%
    dplyr::mutate(call_seg_covered_cr_stable=seg_covered_cr_stable/frac_covered)
  
  sv_cna_fractions_wide = sv_cna_fractions %>% 
    dplyr::select(patient_sv_merged,call,frac_covered,cr_l2fc_50,call_seg_covered_cr_stable) %>% 
    pivot_wider(names_from="call",values_from=c("frac_covered","cr_l2fc_50","call_seg_covered_cr_stable")) 
  
  sv_cna = sv_cna %>% left_join(sv_cna_fractions_wide)
  return(sv_cna)
}

get_sv_cna_flank_overlaps = function(svs,cnas_1mb,svs_df,cnas_df,svs_df_cols=NULL,cnas_df_cols=NULL,sv_id_col="patient_sv_merged", cna_id_col="cna_id"){
  if(is.null(svs_df_cols)) {
    svs_df_overlap_cols = svs_df_overlap_cols[svs_df_overlap_cols %in% names(svs_df)]
    svs_df_anno_cols = svs_df_anno_cols[svs_df_anno_cols %in% names(svs_df)]
    svs_df_cols=c(svs_df_overlap_cols,svs_df_anno_cols)
  } 
  if(is.null(cnas_df_cols)) {
    cnas_df_cols = modeled_seg_cols[modeled_seg_cols %in% names(cnas_df)]
  }
  
  flank_cna_overlaps= get_reciprocal_overlap_pairs_start_end(svs,cnas_1mb,reciprocal_overlap = 0,svtype_matching = FALSE)
  
  if(nrow(flank_cna_overlaps)==0) { return(data.frame()) }
  join_cna_df = cnas_df[,unique(c(cnas_df_cols,"start","end"))] %>% dplyr::rename(cna_start=start,cna_end=end) %>% unique()
  join_svs_df = svs_df[,unique(c(svs_df_cols,"start","end"))] %>% dplyr::rename(sv_start=start,sv_end=end) %>% unique()
  
  flank_cna_overlaps = flank_cna_overlaps %>% dplyr::rename(!!sym(sv_id_col):=set1, !!sym(cna_id_col):=set2, svtype=set1_svtype) %>% 
    left_join(join_cna_df)  %>%
    left_join(join_svs_df) 
  
  if(FALSE) {
    flank_cna_overlaps = flank_cna_overlaps %>% rowwise() %>% mutate(distance = ifelse(sv_breakpoint_orientation=="start",
                                                                                       min(c(abs(sv_start-cna_start),abs(sv_start-cna_end))),
                                                                                       min(c(abs(sv_end-cna_start),abs(sv_end-cna_end))))
    ) %>% as.data.frame() 
  }
  
  
  flank_cna_overlaps = flank_cna_overlaps %>% rowwise() %>% mutate(distance = ifelse(sv_breakpoint_orientation=="start",
                                                                                     ifelse(sv_start>cna_start & sv_start<cna_end,0,  #falls in between
                                                                                            min(c(abs(sv_start-cna_start),abs(sv_start-cna_end)))),
                                                                                     ifelse(sv_end>cna_start & sv_end<cna_end,0,  #falls in between
                                                                                            min(c(abs(sv_end-cna_start),abs(sv_end-cna_end)))))
  ) %>% as.data.frame() 
  
  
  
  flank_cna_overlaps$distance_kbp = flank_cna_overlaps$distance/1e3
  return(flank_cna_overlaps) 
}
annotate_cna_chr_arm = function(cnas,modeled_seg,chr_arms) {
  chr_arm_cna_overlaps = get_reciprocal_overlap_pairs_start_end(cnas,chr_arms,reciprocal_overlap = 0,svtype_matching = F)
  #NEXT: annotate cnas with chr arms and check width compared to full arm
  chr_arm_cna_overlaps = chr_arm_cna_overlaps %>% select(set1,sv_breakpoint_orientation,set2) %>% unique() %>%  
    pivot_wider(names_from="sv_breakpoint_orientation",values_from="set2")
  
  chr_arm_cna_overlaps = chr_arm_cna_overlaps %>% dplyr::rename(cna_id=set1) %>% group_by(cna_id) %>% 
    summarize(chr_arm = ifelse(start!=end, 
                               paste0(start,"-",end),
                               start))
  
  modeled_seg = modeled_seg %>% left_join(chr_arm_cna_overlaps)
  
  return(modeled_seg)
}

get_cna_merged_gainloss = function(cnas,modeled_seg,
                                   cnas_id_col="cna_id",cnas_df_cols = c("cna_id","cr_l2fc_50","maf_50","call","cna_coordinate"),
                                   gap_size=6e3,group_cols="",split_by_cntype=c("gain","loss")) {
  if(is.null(split_by_cntype)) {
    split_by_cntype = cnas$call %>% unique()
  }
  split_by_cntype = split_by_cntype[split_by_cntype %in% unique(cnas$call)]
  ##function often called with already split data so optional to turn of splitting 
  ## todo: check if non overlapping because possible in theory
  cna_merged=GRanges()
  for(cntype in split_by_cntype) {
    cna_gain = cnas[cnas$call==cntype]
    if(length(cna_gain)==0) {
      print(paste0("check cn type ",cntype))
      print(cnas$call %>% unique())
      print(names(cnas) %>% head())
      next()
    }
    cna_gain_merged = GenomicRanges::reduce(cna_gain,min.gapwidth=gap_size) #size of LINE/SINE/SD
    cna_gain_merged$call=cntype
    cna_merged=c(cna_merged,cna_gain_merged)
  }
  
  names(cna_merged)=paste0("cna_merged_",1:length(cna_merged))
 
  cna_merged_df = cna_merged %>% as.data.frame()
  cna_merged_df$cna_merged_id = rownames(cna_merged_df)
  cna_merged_df = cna_merged_df %>% get_gr_coordinate("cna_coordinate")
  
  overlaps_cna_merged = get_reciprocal_overlap_pairs(cnas,cna_merged,reciprocal_overlap = 0,svtype_matching = F)
  
  cna_merged_values = overlaps_cna_merged %>% dplyr::rename(!!sym(cnas_id_col):=set1,cna_merged_id=set2) %>% 
    left_join(modeled_seg[,c(cnas_df_cols)]) %>%
    group_by(across(c("cna_merged_id",group_cols))) %>% 
    summarize(cr_l2fc_50_max = max(cr_l2fc_50,na.rm = T),
              maf_50 = sum(maf_50*overlap_set2_set1,na.rm = T), 
              cr_l2fc_50 = sum(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              cr_sd = sd(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              seg_cnt=n(), 
              loss_cnt=sum(call=="loss"), gain_cnt=sum(call=="gain"), loh_cnt=sum(call=="loh"),
              seg_covered = sum(overlap_set2_set1),
              .groups="drop") 
  cna_merged_df = cna_merged_df %>% left_join(cna_merged_values)
  cna_merged_df = cna_merged_df %>% dplyr::rename(cna_width=width)
  
  
  cna_merged=GRanges(cna_merged_df)
  names(cna_merged) = cna_merged_df$cna_merged_id
  return(cna_merged)
}
annotate_cna_neutral_gaps = function(cna_merged,chromosomes,cnas,modeled_seg) {
  cna_merged_with_neutral=cna_merged
  chr_ends=chromosomes_df[,c("seqnames","end")]
  chr_ends$start=chr_ends$end
  chr_ends = GRanges(chr_ends)
  chr_ends$call="0"
  names(chr_ends)=paste0("chr_ends_",1:length(chr_ends))
  chr_ends$cr_l2fc_50=0
  chr_ends = subsetByOverlaps(chr_ends,cna_merged,invert = T)  # remove already present
  
  cna_merged_with_neutral = c(cna_merged_with_neutral,GRanges(chr_ends))
  gaps_to_add = gaps(cna_merged_with_neutral)
  names(gaps_to_add)=paste0("gap_",1:length(gaps_to_add))
  gaps_df = gaps_to_add %>% as.data.frame()
  gaps_df$gap_id=rownames(gaps_df)
  
  #gaps_to_add$call="0"
  
  overlaps_gaps = get_reciprocal_overlap_pairs(cnas,gaps_to_add,reciprocal_overlap = 0,svtype_matching = F)
  
  gaps_values = overlaps_gaps %>% dplyr::rename(cna_id=set1,gap_id=set2) %>% 
    left_join(modeled_seg[,c(cnas_df_cols)]) %>%
    group_by(gap_id) %>% 
    summarize(cr_l2fc_50_max = max(cr_l2fc_50,na.rm = T),
              maf_50 = sum(maf_50*overlap_set2_set1,na.rm = T), 
              cr_l2fc_50 = sum(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              cr_sd = sd(cr_l2fc_50*overlap_set2_set1,na.rm = T),
              seg_cnt=n(), 
              loss_cnt=sum(call=="loss"), gain_cnt=sum(call=="gain"), loh_cnt=sum(call=="loh"),
              seg_covered = sum(overlap_set2_set1),
              .groups="drop") 
  gaps_df = gaps_df %>% left_join(gaps_values)
  gaps_df = gaps_df %>% dplyr::rename(cna_width=width)
  gaps_df = gaps_df %>% call_cna()
  
  
  gaps_to_add = GRanges(gaps_df %>% dplyr::rename(cna_merged_id = gap_id))
  names(gaps_to_add)=gaps_to_add$cna_merged_id
  cna_merged_with_neutral = c(cna_merged_with_neutral,gaps_to_add)
  
  return(cna_merged_with_neutral)
}

scale_cn_threshold = function(df,threshold=2){
  df = df %>% 
    mutate(cr_l2fc_50= ifelse(cr_l2fc_50>threshold,threshold,cr_l2fc_50)) %>%
    mutate(cr_l2fc_50= ifelse(cr_l2fc_50<(-threshold),-threshold,cr_l2fc_50))
  
  return(df)
}



get_recurrent_cna_count_seg = function(cna_df,group_col="gene_name") {
  recurrent_cna = cna_df %>% filter(gain_cnt>0|loss_cnt>0|loh_cnt) %>% 
    group_by_(group_col) %>%
    summarize(
      patient_cnt = length(unique(patient_id)),
      patient_lst =toString(sort(unique(patient_id))),
      cnl2fc_mean = mean(cna_l2fc_weighted_avg),
      loss_cnt=sum(loss_cnt),gain_cnt=sum(gain_cnt),loh_cnt=sum(loh_cnt)) %>%
    arrange(-patient_cnt)
  return(recurrent_cna)
}

#gene_name	gene_id	gene_type	cr_l2fc_50_max	maf_50	cr_l2fc_50	maf_10	cr_l2fc_10	maf_90	cr_l2fc_90	seg_cnt	loss_cnt	gain_cnt	loh_cnt	seg_covered	overlapping_cna	overlapping_cna_id	patient_id

get_recurrent_gene_cna = function(cna_df,group_cols=c("gene_id","gene_name","call")) {
  recurrent_cna = cna_df %>% 
    group_by(across(all_of(group_cols))) %>%
    summarize(
      patient_cnt = length(unique(patient_id)),
      patient_lst =toString(sort(unique(patient_id))),
      cr_l2fc_mean = mean(cr_l2fc_50),
      loss_cnt=sum(loss_cnt),gain_cnt=sum(gain_cnt),loh_cnt=sum(loh_cnt)) %>%
    arrange(-patient_cnt)
  return(recurrent_cna)
}
