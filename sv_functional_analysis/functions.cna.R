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

if(FALSE){
call_cna = function(cna,cna_cutoff=0.15) {
  convert_gr = F
  if(isS4(cna)) { 
    cna_names = names(cna)
    cna = as.data.frame(cna) 
    convert_gr=T
  }
  cna = cna %>% mutate(call = ifelse(copy_ratio_log2_mean>cna_cutoff,"gain",
                                     ifelse(copy_ratio_log2_mean<(-cna_cutoff),"loss","0")))
  if(convert_gr) {
    cna = GRanges(cna)
    names(cna)=cna_names
  }
  return(cna)
}
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

get_chr_arms = function(chromosome_bands_df) {
  ## note; 1 based coordinates??? 
  
  autosomes= c(paste("chr",1:22,sep=""),"chrX","chrY")
  
  chr_arms = chromosome_bands_df %>% 
    filter(seqnames %in% autosomes) %>% 
    mutate(chr_arm = ifelse(grepl("p",cytoband),paste0(seqnames,"p"),paste0(seqnames,"q"))) %>% 
    mutate(chr_arm = ifelse(grepl("gneg",gieStain)|grepl("gpos",gieStain),chr_arm,paste0(chr_arm,"_",gieStain))) %>% 
    group_by(seqnames,chr_arm) %>% 
    summarize(start = min(start)+1,
              end = max(end))
  return(chr_arms)
}

## from screen rmd

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
