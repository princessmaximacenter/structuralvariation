## Analyse if amplicons connected by SVs
## Last update: 2023-09-05
# input amplicon merged seg and footprints 
# do (low AF) SVs connect amplicons on different chromosomes or multiple peaks within a large footprint.


suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(AnnotationDbi, quietly=TRUE)
  library(VariantAnnotation, quietly=TRUE)
  library(StructuralVariantAnnotation, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  #library(stringdist, quietly=TRUE)
})


wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(wdir,"functions.cna.R"))
source(paste0(wdir,"functions.general.R"))
source(paste0(wdir,"functions.svs.R"))

source(paste0(wdir,"solids.conf"))

# Settings ----

amplicons_footprint_size_min = 50e3 # min size of amplicon to be included in analysis
amplicons_sv_distance_max = 1e3 #to link svs to amplicons
amplicons_distance_apart_min = 5e6 #witin one footprint


## Inputs
run_date=20230701
amplicon_merged_segments_path = paste0(cohort_results_dir,"amplicon_merged_segments.",run_date,".tsv")
amplicon_footprints_path = paste0(cohort_results_dir,"amplicon_footprints.",run_date,".tsv")
merged_svs_classes_path = paste0(cohort_results_dir,"merged_svs_classes.",run_date,".tsv")

## Outputs
amplicon_footprints_annotated_path = paste0(cohort_results_dir,"amplicon_footprints.with_svs.",run_date,".tsv")
svs_connecting_amplicons_path = paste0(cohort_results_dir,"svs_connecting_amplicons.",run_date,".tsv")

merged_svs_cols = c("patient_sv_merged","sv_merged_coordinate","svtype","tumor_af_mean","svlen_mean",
                    "cytoband","partner_cytoband","chrom","partner_chrom","partner_sv_merged")


# Cohort ----
cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

# Resources ----
## Chromosomes ----
map_template_vars=c('${resources_dir}'=resources_dir)
chromosome_bands_path = stri_replace_all_fixed(chromosome_bands_path_template,names(map_template_vars), map_template_vars,vectorize=F)
chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)

chr_arms = get_chr_arms(chromosome_bands_df)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)

chromosome_bands = GRanges(chromosome_bands_df)
names(chromosome_bands)=chromosome_bands$cytoband


#chr_arm_order = chr_arms_df[gtools::mixedorder(chr_arms_df$chr_arm),c("chr_arm")]# %>% flatten_chr()
#chrom_order = chr_arms_df[gtools::mixedorder(chr_arms_df$seqnames),c("seqnames")] %>% unique()# %>% flatten_chr()


# Read  data ----

## Load low tumor AF SVs ----

no_taf_svs_df = load_cohort_svs(cohort,c('${merged_svs_dir}'=merged_svs_dir),svs_path_template=svs_union_anno_somatic_multitool_path_template,flag_add_patient_label = T)

no_taf_svs_df = no_taf_svs_df %>% annotate_sv_cytoband(chromosome_bands)%>% annotate_sv_chr_arm(chr_arms)
no_taf_svs_df = no_taf_svs_df %>% rowwise() %>% mutate(partner_chrom = unlist(strsplit(partner_coordinate,":"))[1][]) %>% as.data.frame()
no_taf_svs_df = annotate_sv_multitool_support(no_taf_svs_df,sv_tumor_af_threshold = 0)

no_taf_merged_svs=make_merged_svs(no_taf_svs_df)
no_taf_merged_svs = no_taf_merged_svs  %>% left_join(cohort[,patient_tumor_id_cols])  

no_taf_merged_svs = no_taf_merged_svs %>% filter(flag_in_filtered_svs) %>% filter(patient_label %in% cohort$patient_label)
no_taf_merged_svs = no_taf_merged_svs %>% filter(chrom %in% autosomes)
resolve_orphan_svs(no_taf_merged_svs,no_taf_svs_df)

## Load SVs > merged svs classes ----
merged_svs = read.table(merged_svs_classes_path,sep="\t",header=T)

merged_svs_gr = GRanges(merged_svs$sv_merged_coordinate)
mcols(merged_svs_gr) = merged_svs %>% select(-start,-end,-seqnames) 
names(merged_svs_gr) = merged_svs_gr$patient_sv_merged


## combine low AF and merged svs because difference in merging can cause slight differences in identifier

all_merged_svs = rbind(no_taf_merged_svs,
                       merged_svs %>% filter(!patient_sv_merged %in% no_taf_merged_svs$patient_sv_merged) %>% select(names(no_taf_merged_svs)))

#annotate here because merged svs overridden by the resolve orphans function
all_merged_svs = all_merged_svs %>% rowwise() %>% dplyr::mutate(chrom_pair=paste0(sort(c(chrom,partner_chrom)),collapse = "-")) %>% 
  left_join(merged_svs %>% select(patient_sv_merged,contains("complex")))  %>% as.data.frame()

all_merged_svs_gr = GRanges(all_merged_svs$sv_merged_coordinate)
mcols(all_merged_svs_gr) = all_merged_svs 
names(all_merged_svs_gr) = all_merged_svs_gr$patient_sv_merged

#always uq
#all_merged_svs %>% group_by(patient_sv_merged) %>% summarize(coord=cnt_str(sv_merged_coordinate)) %>% filter(coord>1) %>% nrow()==0

## Load amplicon segments and footprints ----
amplicon_merged_seg_df = read.table(amplicon_merged_segments_path,sep = "\t",header=T)
amplicon_merged_seg_df = amplicon_merged_seg_df %>% dplyr::rename_with(.cols=c("footprint_id","cna_width","cna_coordinate"),.fn=function(x){paste0("amplicon_",x)})

amplicon_merged_seg = GRanges(amplicon_merged_seg_df$amplicon_cna_coordinate)
mcols(amplicon_merged_seg) = amplicon_merged_seg_df %>% select(-seqnames,-start,-end,-width,-strand)
names(amplicon_merged_seg) = amplicon_merged_seg_df$patient_cna_merged_id

amplicon_footprints = read.table(amplicon_footprints_path,sep = "\t",header=T)
amplicon_footprints$footprint_width=amplicon_footprints$footprint_size

amplicon_footprints = amplicon_footprints %>% dplyr::rename_with(.cols=c("footprint_id","footprint_width","footprint_coordinate"),.fn=function(x){paste0("amplicon_",x)})
amplicon_footprints = amplicon_footprints %>% filter(amplicon_footprint_width>amplicons_footprint_size_min)

# Overlap with SVs ----

map_amplicons_svs = get_reciprocal_overlap_pairs(resize_gr_distance(amplicon_merged_seg,amplicons_sv_distance_max),all_merged_svs_gr,reciprocal_overlap = 0,svtype_matching = F)
map_amplicons_svs = map_amplicons_svs %>% dplyr::rename(patient_cna_merged_id=set1,patient_sv_merged=set2) 

map_amplicons_svs = map_amplicons_svs  %>% select(-from,-to) %>% 
  left_join(amplicon_merged_seg_df[,c("patient_label","patient_cna_merged_id","amplicon_cna_coordinate","cr_l2fc_50","amplicon_footprint_id","amplicon_cna_width")],by="patient_cna_merged_id") %>% 
  left_join(all_merged_svs %>% select(patient_label,patient_sv_merged,sv_merged_coordinate,contains("complex")),by="patient_sv_merged") %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y)

map_amplicons_svs$distance = get_gr_distance(map_amplicons_svs$amplicon_cna_coordinate,map_amplicons_svs$sv_merged_coordinate,flag_as_breakpoints = F) 


#SVs start end

map_amplicons_svs_start_end = get_reciprocal_overlap_pairs_start_end(all_merged_svs_gr,resize_gr_distance(amplicon_merged_seg,amplicons_sv_distance_max),reciprocal_overlap = 0,svtype_matching = F)
map_amplicons_svs_start_end = map_amplicons_svs_start_end %>% dplyr::rename(patient_sv_merged=set1,patient_cna_merged_id=set2) 

map_amplicons_svs_start_end = map_amplicons_svs_start_end  %>% select(-from,-to) %>% 
  left_join(amplicon_merged_seg_df[,c("patient_label","patient_cna_merged_id","amplicon_cna_coordinate","cr_l2fc_50","amplicon_footprint_id","amplicon_cna_width")],by="patient_cna_merged_id") %>% 
  left_join(all_merged_svs %>% select(patient_label,patient_sv_merged,contains("complex")),by="patient_sv_merged") %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y)


svs_amplicons_start_end = map_amplicons_svs_start_end %>% select(c("patient_label","patient_sv_merged","patient_cna_merged_id","sv_breakpoint_orientation")) %>% pivot_wider(id_cols = c("patient_label","patient_sv_merged"),names_from="sv_breakpoint_orientation",values_from = "patient_cna_merged_id")
svs_amplicons_start_end = svs_amplicons_start_end %>% mutate(both_sides=!is.na(start)&!is.na(end))



# Amplicons linked by SVs ----

#amplicons with distance in between on same chrom

amplicon_merged_seg_grl = GenomicRanges::split(amplicon_merged_seg,amplicon_merged_seg$amplicon_footprint_id)
amplicon_merged_seg_distances = get_distance_to_df(amplicon_merged_seg_grl,id_cols = c("patient_label","patient_cna_merged_id","chrom","amplicon_footprint_id"),distance_to_nearest_id_col = c("patient_cna_merged_id"),analyse_sv_bp = F,analysis_type = "consecutive")

amplicon_merged_seg_distances = amplicon_merged_seg_distances %>% left_join(cohort[,c("patient_label","cancer_type")])



#select only svs in amplicon
svs_in_amplicons = all_merged_svs %>% 
  merge(map_amplicons_svs[,c("patient_sv_merged","patient_cna_merged_id")] %>% unique()) %>% 
  left_join( 
    map_amplicons_svs[,c("patient_sv_merged","patient_cna_merged_id")] %>% unique() %>% dplyr::rename(partner_cna_merged_id = patient_cna_merged_id), by=c("partner_sv_merged"="patient_sv_merged"))





## SVs connecting amplicons ----
#different chromosomes, or amplicon of footprints that have large max distance between them

amplicon_distant_pairs = filter(amplicon_merged_seg_distances,upstream_distance>amplicons_distance_apart_min)[,c("patient_cna_merged_id","amplicon_footprint_id","upstream_patient_cna_merged_id","lead_number","upstream_lead_number")] %>% dplyr::rename(end=patient_cna_merged_id,start=upstream_patient_cna_merged_id,end_lead_number=lead_number,start_lead_number=upstream_lead_number) 

amplicon_distant_pairs = rbind(amplicon_distant_pairs,  filter(amplicon_merged_seg_distances,downstream_distance>amplicons_distance_apart_min)[,c("patient_cna_merged_id","amplicon_footprint_id","downstream_patient_cna_merged_id","lead_number","downstream_lead_number")] %>% dplyr::rename(start=patient_cna_merged_id,end=downstream_patient_cna_merged_id,end_lead_number=downstream_lead_number,start_lead_number=lead_number) 
)

amplicon_distant_pairs = amplicon_distant_pairs %>% unique()

#have a look
#amplicon_footprints %>% filter(amplicon_footprint_id %in% amplicon_distant_pairs$amplicon_footprint_id)

#select intra svs connecting the distant amplicons
#intra svs only
intra_svs_connecting_amplicons = svs_amplicons_start_end %>% filter(both_sides & start != end) %>% ungroup() %>% 
  left_join(all_merged_svs) %>% select(merged_svs_cols,contains("complex"),start,end) %>% as.data.frame()

#=> doesnt work because should be "and anything further"
#intra_svs_connecting_amplicons %>% merge(amplicon_distant_pairs)

#add lead number to the intra svs df and look if start lead <= distant pair start and end lead >= distant pair end
intra_svs_connecting_amplicons = intra_svs_connecting_amplicons %>% 
  left_join(amplicon_merged_seg_distances[,c("lead_number","patient_cna_merged_id","amplicon_footprint_id")] %>% dplyr::rename(start=patient_cna_merged_id,start_lead_number=lead_number)) %>%
  left_join(amplicon_merged_seg_distances[,c("lead_number","patient_cna_merged_id")] %>% dplyr::rename(end=patient_cna_merged_id,end_lead_number=lead_number))

intra_svs_connecting_amplicons = intra_svs_connecting_amplicons %>% 
  left_join(amplicon_distant_pairs %>% select(amplicon_footprint_id,start_lead_number,end_lead_number) %>% dplyr::rename(dp_start_lead_number=start_lead_number,dp_end_lead_number=end_lead_number)) 

intra_svs_connecting_amplicons = intra_svs_connecting_amplicons %>%
  filter(start_lead_number <= dp_start_lead_number & end_lead_number >= dp_end_lead_number)


#interchromosomal 
inter_chrom_svs_connecting_amplicons =  svs_in_amplicons %>% filter(partner_sv_merged %in% svs_in_amplicons$patient_sv_merged) %>%
  dplyr::rename(start=patient_cna_merged_id,end=partner_cna_merged_id) %>% select(merged_svs_cols,contains("complex"),start,end)

#combine intra and inter
svs_connecting_amplicons = rbind(intra_svs_connecting_amplicons %>% select(names(inter_chrom_svs_connecting_amplicons)),inter_chrom_svs_connecting_amplicons)

svs_connecting_amplicons = svs_connecting_amplicons %>% 
  left_join(amplicon_merged_seg_df[,c("patient_label","patient_cna_merged_id","amplicon_footprint_id")] %>% dplyr::rename(start_footprint_id=amplicon_footprint_id) ,by=c("start"="patient_cna_merged_id")) %>%
  left_join(amplicon_merged_seg_df[,c("patient_cna_merged_id","amplicon_footprint_id")] %>% dplyr::rename(end_footprint_id=amplicon_footprint_id) ,by=c("end"="patient_cna_merged_id"))


# Integrate findings ----

#multiple peaks within 1 footprint connected or not? 
#todo: add complex sv

#overview of multiple amplicons with max consecutive distance 
#these are cases to check if SVs connect them.

amplicon_footprints_max_consecutive_distance = amplicon_merged_seg_distances %>% filter(!is.na(upstream_distance)) %>% 
  group_by(amplicon_footprint_id) %>% summarize(max_consecutive=max(upstream_distance,na.rm = T),max_consecutive_mbp=max_consecutive/1e6) 

#add which ones
amplicon_footprints_max_consecutive_distance = amplicon_footprints_max_consecutive_distance %>% 
  left_join(amplicon_merged_seg_distances[,c("amplicon_footprint_id","upstream_distance","upstream_patient_cna_merged_id","patient_cna_merged_id")] %>% unique(),by=c("amplicon_footprint_id","max_consecutive"="upstream_distance"))

amplicon_footprints_max_consecutive_distance = amplicon_footprints_max_consecutive_distance %>% 
  dplyr::mutate(connected_by_svs=(amplicon_footprint_id %in% svs_connecting_amplicons$start_footprint_id | amplicon_footprint_id %in% svs_connecting_amplicons$end_footprint_id))



# Both breakpoints in amplicon? ----

distant_map_amplicons_svs_start_end = get_reciprocal_overlap_pairs_start_end(merged_svs_gr,resize_gr_distance(amplicon_merged_seg,100e6),reciprocal_overlap = 0,svtype_matching = F)
distant_map_amplicons_svs_start_end = distant_map_amplicons_svs_start_end %>% dplyr::rename(patient_sv_merged=set1,patient_cna_merged_id=set2) 

distant_map_amplicons_svs_start_end = distant_map_amplicons_svs_start_end  %>% select(-from,-to) %>% 
  left_join(amplicon_merged_seg_df[,c("patient_label","patient_cna_merged_id","amplicon_cna_coordinate","cr_l2fc_50","amplicon_footprint_id","amplicon_cna_width","start","end")],by="patient_cna_merged_id") %>% 
  merge(merged_svs %>% select(patient_sv_merged,contains("complex"),patient_label,chrom,start,end),by="patient_sv_merged") %>%
  filter(patient_label.x==patient_label.y) %>% dplyr::rename(patient_label=patient_label.x) %>% select(-patient_label.y)

#parse one as bp and other as range
distant_map_amplicons_svs_start_end = distant_map_amplicons_svs_start_end %>% mutate(sv_bp_coord = paste0(chrom,":",
                                                                                                          ifelse(sv_breakpoint_orientation=="start",start.y,end.y)))

distant_map_amplicons_svs_start_end$distance = get_gr_distance(distant_map_amplicons_svs_start_end$amplicon_cna_coordinate,distant_map_amplicons_svs_start_end$sv_bp_coord,flag_as_breakpoints = F) 


distant_map_amplicons_svs_start_end$patient_sv_merged_bp = paste0(distant_map_amplicons_svs_start_end$patient_sv_merged,"_",distant_map_amplicons_svs_start_end$sv_breakpoint_orientation)

distant_map_amplicons_svs_start_end = distant_map_amplicons_svs_start_end[distant_map_amplicons_svs_start_end$distance == ave(distant_map_amplicons_svs_start_end$distance,distant_map_amplicons_svs_start_end$patient_sv_merged_bp,FUN=min),]

distant_map_amplicons_svs_start_end_wide = distant_map_amplicons_svs_start_end  %>% select(c("patient_label","patient_sv_merged","patient_cna_merged_id","sv_breakpoint_orientation","distance")) %>% pivot_wider(id_cols = c("patient_label","patient_sv_merged"),names_from="sv_breakpoint_orientation",values_from = c("patient_cna_merged_id","distance"),names_sep = "_")

#one bp inside amplicon and other far away could indicate HSR
#count inside or outside of amplicon relative to amplicons_sv_distance_max

amplicon_svs_position = merged_svs %>% filter(complex_sv_class=="amplicon") %>% 
  select(patient_label,patient_sv_merged, complex_sv_id) %>% 
  left_join(distant_map_amplicons_svs_start_end_wide) %>% 
  mutate(position_relative_to_amplicon = ifelse( (distance_start < amplicons_sv_distance_max & distance_end<amplicons_sv_distance_max), "both_bp_inside",
                                                 ifelse( (distance_start<amplicons_sv_distance_max & (is.na(distance_end) | distance_end>amplicons_sv_distance_max)) |
                                                           distance_end<amplicons_sv_distance_max & (is.na(distance_start) | distance_start>amplicons_sv_distance_max) , "one_bp_inside","outside")),
         position_relative_to_amplicon = ifelse(is.na(position_relative_to_amplicon) & is.na(distance_start) & is.na(distance_end),"outside",position_relative_to_amplicon)) %>%
  group_by(patient_label,complex_sv_id,position_relative_to_amplicon) %>% summarize(cnt=cnt_str(patient_sv_merged)) 

amplicon_svs_position_wide = amplicon_svs_position %>% pivot_wider(id_cols=c("complex_sv_id","patient_label"),names_from = "position_relative_to_amplicon",values_from = "cnt")
amplicon_svs_position_wide[is.na(amplicon_svs_position_wide)]=0




## Amplicons per patient - complex SV yes/no. - multiple Y/N ----

## add max consecutive distance
## add if svs connect amplicons
#the ones in svs_connecting_amplicons are also from the sufficient distance list so can use that as well to infer
svs_connecting_amplicons_summary = svs_connecting_amplicons %>% group_by(start_footprint_id) %>% 
  summarize(connected_footprints_lst=lst_str(c(start_footprint_id,end_footprint_id))) %>% 
  dplyr::mutate(connected_footprints_lst = ifelse(connected_footprints_lst==start_footprint_id,"within",connected_footprints_lst))

amplicon_footprints = amplicon_footprints %>%  
  left_join(amplicon_footprints_max_consecutive_distance %>% select(amplicon_footprint_id,max_consecutive)) %>%
  left_join(svs_connecting_amplicons_summary,by=c("amplicon_footprint_id"="start_footprint_id")) %>%
  dplyr::mutate(multiple_connected_by_svs=!is.na(connected_footprints_lst))

amplicon_footprints = amplicon_footprints %>% dplyr::rename(max_consecutive_distance_within_amplicon = max_consecutive) 
amplicon_footprints[is.na(amplicon_footprints$max_consecutive_distance_within_amplicon),c("max_consecutive_distance_within_amplicon")] =0

#check if all amplicon type
filter(map_amplicons_svs,flag_is_complex&complex_sv_class!="amplicon")$amplicon_footprint_id %>% unique() %>% length() ==0

## add total number of amplicons > possibility of connecting

amplicon_footprints_per_patient = amplicon_footprints %>% group_by(patient_label) %>% summarize(amplicons_lst=lst_str(amplicon_footprint_id),amplicons_cnt=cnt_str(amplicon_footprint_id))

amplicon_footprints = amplicon_footprints %>% left_join(amplicon_footprints_per_patient)
amplicon_footprints = amplicon_footprints %>% mutate(multiple_amplicon_regions_to_connect = amplicons_cnt>1 | max_consecutive_distance_within_amplicon>amplicons_distance_apart_min)


##number of amplified bases
amplified_bases = amplicon_merged_seg_df %>% group_by(patient_label,chrom,amplicon_footprint_id) %>% summarize(amplified_bases=sum(amplicon_cna_width))
amplicon_footprints = amplicon_footprints %>% left_join(amplified_bases)



## add complex svs 
complex_svs_amplicons_summary = map_amplicons_svs %>% filter(flag_is_complex) %>% 
  left_join(amplicon_svs_position_wide) %>%
  group_by(patient_label,amplicon_footprint_id) %>% 
  summarize(complex_svs_lst=lst_str(complex_sv_id),
            complex_svs_bp_position=lst_str(paste0("#svs with both bp inside amplicon: ",both_bp_inside,"; one bp inside: ",one_bp_inside,"; both bp outside: ",outside)))


amplicon_footprints = amplicon_footprints %>% 
  left_join(complex_svs_amplicons_summary) %>%
  dplyr::mutate(overlap_complex_sv=!is.na(complex_svs_lst))


# Export ----

write.table(amplicon_footprints,amplicon_footprints_annotated_path,sep = "\t",col.names = T,row.names = F)
write.table(svs_connecting_amplicons,svs_connecting_amplicons_path,sep = "\t",col.names = T,row.names = F)



