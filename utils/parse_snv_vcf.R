## Parse snv table
## Last update: 2021-12-19
### Changelog
## Refactor for HPC: run for single patient / input file
#should be able to run for vcf input

if(FALSE) {
  
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.conf")
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/sv_functional_analysis/default.docker.conf")
  
  #analysis_type="somatic"
  analysis_type="germline"
  
  run_tumor_id="PMABM000HAS"
  run_normal_id="PMABM000HBB"
 
  
  # run_tumor_id="" #required unless germline
  #  run_normal_id="" #required
  #  run_patient_id=""   #add patient id but do not rely on it being there?.. not sure if sensible
  #  snv_vcf_path = "" #sys glob based on biomaterial id if not provided
  
  source("/hpc/pmc_gen/ivanbelzen/structuralvariation/utils/parse_snv_vcf.R")
  
  
  snv_somatic_panel_path_template="${resources_dir}hg38_diagnostic_somatic_2.4.bed.gz"
}

snv_data_dir="/hpc/pmc_gen/ivanbelzen/structuralvariation/tmp/snv/"
  
blacklist_path="/hpc/pmc_gen/ivanbelzen/wdl_pipeline/fusion_pilot/hg38-blacklist.v2.bed"

#sys glob with biomaterial id if no path provided 
if(!exists("run_normal_id")&(analysis_type!="somatic" | exists("run_tumor_id"))) {
  print("No biomaterial ids provided, exiting")
  #quit()
}

#patient basename tumor_id + _ + normal_id for somatic and normal_id only for germline
if(analysis_type=="somatic"){ 
  patient_basename = paste0(run_tumor_id,"_",run_normal_id)
} else {
  patient_basename = paste0(run_normal_id)
}

if(!exists("snv_vcf_path")) {
  snv_vcf_path = paste0(snv_data_dir,patient_basename,snv_vcf_file_ext)
}

suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(dplyr)
  #library(AnnotationHub)
  
})
source(paste0(script_dir,"functions.general.R"))
source(paste0(script_dir,"functions.svs.R"))
source(paste0(utils_script_dir,"functions.snv.R"))


#patient basename ${tumor_id}_${normal_id} for somatic and ${normal_id} only for germline
map_template_vars=c('${resources_dir}'=resources_dir,'${utils_output_dir}'=utils_output_dir,'${patient_basename}'=patient_basename,'${analysis_type}'=analysis_type)

snv_output_path = stri_replace_all_fixed(snv_output_path_template,names(map_template_vars), map_template_vars,vectorize=F)
snv_impact_output_path = stri_replace_all_fixed(snv_impact_output_path_template,names(map_template_vars), map_template_vars,vectorize=F)

if(analysis_type=="somatic"){ 
  gene_panel_file = stri_replace_all_fixed(snv_somatic_panel_path_template,names(map_template_vars), map_template_vars,vectorize=F)
} else {
  gene_panel_file = stri_replace_all_fixed(snv_germline_panel_path_template,names(map_template_vars), map_template_vars,vectorize=F)
}  



#settings
#snv_overview_cols=c("seqnames", "start", "end", "width", "strand", "paramRangeID", "REF", "ALT", "QUAL", "FILTER", 
#                    "snv_id", "Allele", "VARIANT_CLASS", "Gene", "Feature", "SYMBOL", "CCDS", 
#                    "STRAND", "IMPACT", "Consequence", "clinvar_rs", "Diag_Somatic_Gene","diagnostic_panel",
#                    "tumor_AF","normal_AF","tumor_AD","normal_AD","tumor_DP","normal_DP")

#disable: just get all cols!

## Note: anno hub is disabled for now on HPC
#hub = AnnotationHub()
#ahDb = query(hub, c("EnsDb", "apiens", "92"))
#ahEdb = ahDb[[1]]

## Prepare VCF file

snv_file=Sys.glob(snv_vcf_path)
if(length(snv_file)!=1){
  print(paste0("ERROR: missing file: ",snv_vcf_path))
  quit()
}


## Load somatic panel
gene_panel = rtracklayer::import(gene_panel_file, format="bed")
#note that score column is a placeholder
gene_panel_snv_lst = sort(unique(gene_panel$name))

##Input: SNV file location 
#pre filter on PASS variants only 
if(length(Sys.glob(paste0(snv_file,".tbi")))==1){
  prefilterPASS = function(x) { grepl("PASS", x, fixed=TRUE) }
  pre = FilterRules(list(prefilterPASS = prefilterPASS))
  vcf_filtered = filterVcf(snv_file, "hg38", tempfile(), prefilters=pre)
  snv_vcf = readVcf(vcf_filtered, "hg38")
} else {
  snv_vcf = readVcf(snv_file, "hg38")
  snv_id_pass = names(rowRanges(snv_vcf)[rowRanges(snv_vcf)$FILTER=="PASS"])
  snv_vcf = snv_vcf[snv_id_pass]
}  


    ## Get information from VEP columns, like ensembl gene id, gene name (symbol), predicted effect etc
    # either make you own list of VEP columns or use this snippet to get from the VCF header
    # VEP data is in stored in "CSQ" column of VCF and separated with pipes (|)
    # there are some quirky things like an empty column and a duplicated one according to VCF header 
    info_vcf = as.data.frame(info(header(snv_vcf)))
    vep_column_charstring = str_split(info_vcf["CSQ",c("Description")],"Format: ")[[1]][2]
    vep_columns = str_split(vep_column_charstring,"\\|")[[1]]
    vep_columns[vep_columns==""]="empty" #prevent error
    vep_columns[duplicated(vep_columns)]=paste0(vep_columns[duplicated(vep_columns)],".v2")
    ##if you dont subset, you have to change the duplicated names
    
    snv_vep_df = as.data.frame(info(snv_vcf)) %>% mutate(snv_id = rownames(.)) 
    
    #update 2022-03-09 also get other cols
    #%>% dplyr::select(snv_id,CSQ)
    
    #testing purposes only 
    #snv_vep_df=head(snv_vep_df,n=1000)
    #snv_vep_df = snv_vep_df %>% separate(col="CSQ",into=vep_columns,sep="\\|",extra = "drop")
    
    #update 2023-08-21 expand list to get all predictions and solve missing mutations
    
    snv_vep_df = snv_vep_df %>% unnest_longer(col=c("CSQ")) %>% as.data.frame() #%>% dplyr::select(snv_id,CSQ)
    snv_vep_df = snv_vep_df %>% separate(col="CSQ",into=vep_columns,sep="\\|",extra = "drop") %>% unique()
    
    ## Add from geno() the AF, AD, DP
    snv_vep_df = snv_vep_df[,!names(snv_vep_df) %in% c("AF","AD","DP")]
    
if(analysis_type=="somatic"){ 
  snv_geno_df= parse_snv_geno_somatic(snv_vcf,run_tumor_id,run_normal_id)
} else {
  snv_geno_df= parse_snv_geno_germline(snv_vcf)
}

    snv_vep_df = snv_vep_df %>%  left_join(snv_geno_df,by="snv_id")

    snv_gene_panel = subsetByOverlaps(rowRanges(snv_vcf),gene_panel)
    
    # access snv ids as names(snv_gene_panel)
   # snv_id_selection = unique(c(names(snv_gene_panel), snv_vep_df[snv_vep_df$SYMBOL %in% gene_panel_snv_lst,c("snv_id")]))
    
    #note: gene panel should be only based on intersect not gene names
    snv_vep_df = snv_vep_df %>% mutate(diagnostic_panel = snv_id %in% names(snv_gene_panel))
    
    ## FILTER TO PASS 

    if(FALSE) { 
    ## Summarize per gene fist
    if(analysis_type=="somatic") {
      
  snv_gene_summary = snv_vep_df %>% dplyr::filter(Gene!="" & tumor_AF>0.05) %>% group_by(Gene,SYMBOL) %>% 
    summarize(vep_impact = toString(IMPACT), vep_consequence = toString(Consequence),
                snvs = toString(snv_id), snv_cnt = n(), diagnostic_panel = toString(diagnostic_panel),
                tumor_af_lst = toString(tumor_AF),normal_af_lst = toString(normal_AF))
    } else {
      
      snv_gene_summary = snv_vep_df %>% dplyr::filter(Gene!="" & normal_AF>0.05) %>% group_by(Gene,SYMBOL) %>% 
        summarize(vep_impact = toString(IMPACT), vep_consequence = toString(Consequence),
                  snvs = toString(snv_id), snv_cnt = n(), #diagnostic_panel = toString(diagnostic_panel),
                  normal_af_lst = toString(normal_AF)) %>%
        as.data.frame()
      
    }
   write.table(snv_gene_summary,paste0(snv_summary_dir,snv_gene_summary_outfile,patient$patient_identifier,".tsv"),quote = T,sep = "\t",row.names=FALSE)
    }
    
    ## Filter on gene panel genes
    #rigorous filter cant go back
    #snv_vcf = snv_vcf[snv_id_selection]
    # NOTE: 2021-12-19 dont want do this! do want to filter on pass / by default and tumor AF > 0.05 but can also be hindsight
    
    if(analysis_type=="somatic") {
      snv_id_selection=snv_vep_df[snv_vep_df$tumor_AF>0.05,c("snv_id")]
    } else {
      snv_id_selection=snv_vep_df[snv_vep_df$normal_AF>0.05,c("snv_id")]
    }
    
    snv_id_selection = unique(snv_id_selection)
    
    snv_output_df = as.data.frame(rowRanges(snv_vcf[snv_id_selection]))
    
    alt = as.character(unlist(snv_vcf$ALT))
    if(length(alt)==nrow(snv_output_df)){
      snv_output_df$ALT = alt
    } else {
      snv_output_df$ALT = NA
    }
      
    snv_output_df= snv_output_df %>% mutate(snv_id = rownames(.)) %>% left_join(snv_vep_df,by="snv_id")
    
    ##2022-03-09
    ##disable subset by overview cols -> output all unless
    #snv_overview_cols=snv_overview_cols[snv_overview_cols %in% names(snv_output_df)]
    
    excl_cols_lst=c("Location","Diag_Somatic_Gene","Diag_Somatic_Gene_..fork")
    snv_overview_cols=names(snv_output_df)[!names(snv_output_df) %in% excl_cols_lst]
    snv_output_df = snv_output_df[,snv_overview_cols]
    
    
   if(FALSE) {
    ## TODO: annotation hub does not yet work on HPC
   
    #map to cds only for those in somatic panel 
  cds_selection = snv_output_df[snv_output_df$IMPACT %in% c("HIGH","MODERATE") | snv_output_df$snv_id %in% names(snv_gene_panel),c("snv_id")]
    
   snv_cds_mapping = data.frame()
  for(snv_name in cds_selection){
      snv = rowRanges(snv_vcf)[snv_name]
      
      genome(snv)="GRCh38"
      seqlevelsStyle(snv)="Ensembl"
      
      tx_locations = genomeToTranscript(snv,ahEdb)
      tx_locations = unlist(tx_locations)
      if(tx_locations@start==-1){next()}
      tx_locations = transcriptToCds(tx_locations,ahEdb)
      
      properties = snv_output_df[snv_output_df$snv_id==snv_name,]

      snv_location = tx_locations[mcols(tx_locations)$tx_id == properties$Feature,]
      snv_location = as.data.frame(snv_location)
      if(nrow(snv_location)>0){
      gene_prop = genes(x = ahEdb,filter = ~gene_id == properties$Gene)
      snv_location$snv_id = snv_name
      snv_location$seq_strand = as.character(gene_prop@strand)
      if(snv_location$seq_strand=="-"){
        #note: keep the toString because it handles the DNAStringSet object well and as.character does not
        snv_location$ref = toString(complement(mcols(snv)$REF))
        snv_location$alt = toString(complement(unlist(mcols(snv)$ALT)))
      } else {
        snv_location$ref = toString(mcols(snv)$REF)
        snv_location$alt = toString(unlist(mcols(snv)$ALT))
      }
      snv_location$gene_name = properties$SYMBOL
      snv_location$ensembl_id = properties$Gene
      
      snv_location$cds_mapping = paste0(snv_location$gene_name,"_c.",snv_location$start,snv_location$ref,">",snv_location$alt)
      
      
      snv_cds_mapping = rbind(snv_cds_mapping,snv_location)
      }
    }

  if(nrow(snv_cds_mapping)>0){
      write.table(snv_cds_mapping,paste0(snv_summary_dir,snv_cds_mapping_outfile,patient$patient_identifier,".tsv"),quote = T,sep = "\t",row.names=FALSE)
      
      snv_output_df = snv_output_df %>% left_join(snv_cds_mapping[,c("snv_id","cds_mapping")],by="snv_id")
    } else {
    snv_output_df$cds_mapping=NA
    }
   }
    
  ## get rid of sometimes linebreak in between?
    
    if("Diag_Somatic_Gene" %in% names(snv_output_df)) {
    snv_output_df$Diag_Somatic_Gene = gsub("[\r\n]", "",snv_output_df$Diag_Somatic_Gene)
    }


    if("Diag_Germline_Gene" %in% names(snv_output_df)) {
    snv_output_df$Diag_Germline_Gene = gsub("[\r\n]", "",snv_output_df$Diag_Germline_Gene)
    }
    
    if(exists("run_patient_id")) { snv_output_df$patient_id=run_patient_id }

## flag overlap blacklist
    if(length(Sys.glob(blacklist_path))==1){
    blacklist = rtracklayer::import.bed(blacklist_path)
    
    snv_gr = GRanges(snv_output_df[,c("start","end","seqnames","snv_id")])
    names(snv_gr) = snv_gr$snv_id
    snv_in_blacklist = subsetByOverlaps(snv_gr,blacklist)
    
    snv_output_df = snv_output_df %>% mutate(dac_blacklist= snv_id %in% snv_in_blacklist$snv_id)
    }

# writing comment with VEP and diagnostic panel version
vep_version = as.data.frame(meta(header(snv_vcf))$VEP)$Value

if(analysis_type=="somatic"){
  
con <- file(snv_output_path, open="wt")
writeLines(paste0("#",vep_version), con)
writeLines(paste0("#",gene_panel_file), con)
write.table( snv_output_df , con,quote = T,sep = "\t",row.names=FALSE)
close(con)

}
#for germline only save high/moderate

con <- file(snv_impact_output_path, open="wt")
writeLines(paste0("#",vep_version), con)
writeLines(paste0("#",gene_panel_file), con)
write.table( snv_output_df %>% filter(IMPACT %in% c("HIGH","MODERATE")) , con, quote = T,sep = "\t",row.names=FALSE)
close(con)


    
