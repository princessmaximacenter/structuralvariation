# Oncoprint
# Last update 2023-09-05 
# 
# todo refactor clinrel_mutations_incl_prognostic => driver_gene_alterations
# todo paths 
# todo description

## Settings ----
library(GenomicRanges)
library(tidyverse, quietly=TRUE)
library(stringr, quietly=TRUE)
library(testthat, quietly=TRUE)
library(stringdist, quietly=TRUE)
#library(VennDetail)
library(ggplot2)
library(RColorBrewer)
library(DiagrammeR)
library(ggrepel)
library(GGally)
library(pheatmap)
library(dplyr)
library(stringi)
library(lemon)

wdir="~/PycharmProjects/structuralvariation/sv_functional_analysis/"
source(paste0(wdir,"default.conf"))
source(paste0(wdir,"functions.cna.R"))
source(paste0(wdir,"functions.general.R"))
source(paste0(wdir,"functions.svs.R"))

source(paste0(wdir,"solids.conf"))

plot_dir=cohort_results_dir

run_date=20230701

#input
curated_complex_sv_cn_change_path = paste0(cohort_results_dir,"map_complex_sv_cn_change.",run_date,".curated.tsv")
map_clinrel_chromalt_path = paste0(resources_dir,"map_clinrel_chromalt.txt")
driver_gene_alterations_path = paste0(cohort_results_dir,"clinrel_mutations_incl_prognostic.",run_date,".tsv")
ctx_chrom_cn_properties_path = paste0(cohort_results_dir,"ctx_chrom_cn_properties.",run_date,".tsv")
chromosome_cn_calls_path = paste0(cohort_results_dir,"chromosome_level_cn_calls.",run_date,".tsv")

#output
oncoprint_data_path = paste0(cohort_results_dir,"oncoprint_data.",run_date,".tsv")



# Cohort ----
#Load cohort and annotate 

cohort = read.table(patient_table_path,sep = "\t", header=T) %>% arrange(patient_id)
patient_id_to_labels = cohort[,c("patient_id","patient_label")]
patient_tumor_id_cols = c("patient_id","patient_label","tumor_id","cancer_type")

cohort %>% filter(!cancer_type  %in% names(annotation_colors$cancer_type)) %>% select(cancer_type) %>% unique()



# Read files ----
## Load patient vs gene overview ----


patient_gene_overview_detailed_cancer_path = paste0(cohort_results_dir,"gene_level_overview.somatic.detailed.cancer.tsv")
patient_gene_overview_all_cancer = read.table(patient_gene_overview_detailed_cancer_path,sep="\t",header=T)
patient_gene_overview_all_cancer = patient_gene_overview_all_cancer %>% mutate(alteration_simple=trimws(alteration_simple)) 
patient_gene_overview_all_cancer = patient_gene_overview_all_cancer %>% dplyr::mutate(patient_gene = paste0(patient_label,"_",gene_name)) 



## load curated overviews clinrel mutations and CN changes ----

clinrel_mutations_incl_prognostic = read.table(driver_gene_alterations_path, sep="\t",header=T)
map_complex_sv_cn_change = read.table(curated_complex_sv_cn_change_path,sep="\t",header = T)
map_clinrel_chromalt = read.table(map_clinrel_chromalt_path,sep="\t",header=T)

map_clinrel_chromalt = map_clinrel_chromalt %>% mutate(chr_arm_cn_change=paste(chr_arm,cn_change,sep = "_"))



## load unbalanced translocations and aneuploidy ----
ctx_chrom_cn_properties =read.table(ctx_chrom_cn_properties_path,header=T,sep="\t")
unbalanced_ctx = ctx_chrom_cn_properties %>% filter(chrom_cn_unbalanced)
unbalanced_ctx = unbalanced_ctx %>% left_join(cohort[,c("patient_label","cancer_type")])

#klopt niet helemaal, kunnen er meer zijn!
#unbalanced_ctx$chrom_cn_judgement
unbalanced_ctx = unbalanced_ctx %>%
 mutate(chr_arm_cn_change=paste(chr_arm,chrom_cn_judgement,sep = "_"),
         patient_chr_alt=paste(patient_label,chr_arm,chrom_cn_judgement,sep = "_"))

clinrel_unbalanced_ctx = unbalanced_ctx %>% merge(map_clinrel_chromalt) %>% filter(chrom_cn_judgement==cn_change)
unbalanced_ctx = unbalanced_ctx %>% mutate(call=chrom_cn_judgement,clinrel=patient_chr_alt %in% clinrel_unbalanced_ctx$patient_chr_alt)

#remove already in complex
#unbalanced_ctx %>% filter(patient_chr_alt %in% map_complex_sv_cn_change$patient_chr_alt)
unbalanced_ctx = unbalanced_ctx %>% filter(!patient_chr_alt %in% map_complex_sv_cn_change$patient_chr_alt)


# aneuploidy

#fullchrom seqnames as 2 rows 
#removed unbalanced ctx
chromosome_cn_calls = read.table(chromosome_cn_calls_path,sep="\t",header = T)

numerical_cn_calls = chromosome_cn_calls %>%  
  filter(seqnames %in% autosomes) %>% #includes X
  filter(!flag_overlap_unb_ctx)  %>% filter(call!="0" & selected)

numerical_cn_chr_arms = numerical_cn_calls %>% filter(region=="chr_arm") %>% select(patient_label,seqnames,chr_arm,call,cancer_type,region) 
numerical_cn_chr_arms = rbind(numerical_cn_chr_arms,
                              numerical_cn_calls %>% filter(region=="chrom") %>% select(patient_label,seqnames,call,cancer_type,region) %>%
                                mutate(chr_arm=paste0(seqnames,"p")))
numerical_cn_chr_arms = rbind(numerical_cn_chr_arms,
      numerical_cn_calls %>% filter(region=="chrom") %>% select(patient_label,seqnames,call,cancer_type,region) %>% mutate(chr_arm=paste0(seqnames,"q")))

numerical_cn_chr_arms = numerical_cn_chr_arms %>%
 mutate(chr_arm_cn_change=paste(chr_arm,call,sep = "_"),
         patient_chr_alt=paste(patient_label,chr_arm,call,sep = "_"))


clinrel_numerical_cn_chr_arms = numerical_cn_chr_arms %>% merge(map_clinrel_chromalt) %>% filter(call==cn_change)

numerical_cn_chr_arms = numerical_cn_chr_arms %>% mutate(clinrel=patient_chr_alt %in% clinrel_numerical_cn_chr_arms$patient_chr_alt)

#remove already in complex
numerical_cn_chr_arms %>% filter(patient_chr_alt %in% unbalanced_ctx$patient_chr_alt)  #should be 0
numerical_cn_chr_arms %>% filter(patient_chr_alt %in% map_complex_sv_cn_change$patient_chr_alt)

numerical_cn_chr_arms = numerical_cn_chr_arms %>% filter(!patient_chr_alt %in% map_complex_sv_cn_change$patient_chr_alt)
## should use GR overlaps instead of this?
#example both:
#       M959AAB     chr2   chr2p gain              neuroblastoma   chrom        chr2p_gain  M959AAB_chr2p_gain    TRUE
#       M962AAB     chr2   chr2p gain              neuroblastoma   chrom        chr2p_gain  M962AAB_chr2p_gain    TRUE
#M911AAA chr6p loss due to complex sv indeed 




# Figures ----
## colors and patient order ----

extra_annotation_labels = cohort[,c("patient_label","patient_id","cancer_type")]
extra_annotation_labels$gene_name="annotation"

label_order_complex_driver = read.table(paste0(cohort_results_dir,"label_order_complex_driver.tmp")) %>% flatten_chr()
label_order_complex_effect= read.table(paste0(cohort_results_dir,"label_order_complex_effect.tmp")) %>% flatten_chr()

cancer_type.labs = c("EWS","NBL","fpRMS","fnRMS","WT","HBL")
names(cancer_type.labs) = c("ewing_sarcoma","neuroblastoma","alveolar_rhabdomyosarcoma","embryonal_rhabdomyosarcoma","nephroblastoma","hepatoblastoma")

annotation_detail_colors = c(`complex`="red3",`aneuploidy`="mediumpurple3",annotation_colors$alteration_simple)
annotation_detail_colors_incl_aneuploidy = c(`complex`="red3",`aneuploidy`="#e5c494",`cna` = "#e5c494",`snv` = "mediumpurple2", `sv` = "#a6d854")

plot_version = 20230828



## Driver genes (clinrel_mutations_incl_prognostic) ----

# todo refactor clinrel_mutations_incl_prognostic => driver_gene_alterations
gene_order_complex =  clinrel_mutations_incl_prognostic %>% filter(complex_sv_alteration) %>% get_recurrent_gene_centric(flag_patient_label = T) %>% arrange(-patient_cnt) %>% select(gene_name) %>% flatten_chr()

gene_order =  patient_gene_overview_all_cancer %>% filter(
  flag_affected_by_snv | flag_affected_by_sv | flag_cn_hom_loss | flag_cn_amplification |  flag_affected_by_sv_1mb_cn_change) %>% filter(gene_name %in% clinrel_mutations_incl_prognostic$gene_name) %>% get_recurrent_gene_centric(flag_patient_label = T) %>% arrange(-patient_cnt) %>% select(gene_name) %>% flatten_chr()

gene_order= c(gene_order,"annotation") %>% unique()

gene_order_incl_complex = c(gene_order_complex,gene_order,"annotation") %>% unique()

plot_display_gene_overview = patient_gene_overview_all_cancer %>% filter(
  flag_affected_by_snv | flag_affected_by_sv | flag_cn_hom_loss | flag_cn_amplification |  flag_affected_by_sv_1mb_cn_change) %>% 
  filter(gene_name %in% clinrel_mutations_incl_prognostic$gene_name)

#add manual from clinrel mutations incl prognostic such as germline tp53
manual_add = clinrel_mutations_incl_prognostic %>% filter(!patient_gene %in% plot_display_gene_overview$patient_gene)
plot_display_gene_overview = rbind(plot_display_gene_overview,manual_add[,names(plot_display_gene_overview)])

plot_display_gene_overview = rbind_no_colmatch(plot_display_gene_overview,extra_annotation_labels)

plot_display_gene_overview = plot_display_gene_overview %>%
  left_join(clinrel_mutations_incl_prognostic %>%
              select(patient_label,gene_name,driver_gene,complex_sv_alteration) %>% mutate(clinrel_or_prognostic=T) %>% unique(),
            by=c("patient_label","gene_name"))


plot_display_gene_overview[is.na(plot_display_gene_overview$driver_gene),c("driver_gene")]=F
plot_display_gene_overview[is.na(plot_display_gene_overview$clinrel_or_prognostic),c("clinrel_or_prognostic")]=F

plot_display_gene_overview = plot_display_gene_overview %>% 
  dplyr::mutate(alteration_display = 
 ifelse(!is.na(complex_sv_alteration)&complex_sv_alteration,"complex",
 ifelse(grepl("sv",alteration_simple),"sv",
 ifelse(alteration_simple %in% c("amp","hom_loss"), "cna", alteration_simple))))


#gene_order
plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(gene_order_incl_complex))


#order_patient_svs_cnt
#label_order_complex_driver
plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_effect))


display_values = plot_display_gene_overview$alteration_display %>% unique()
display_values[!display_values %in% names(annotation_detail_colors) ]


plot_display_gene_overview$cancer_type = factor(plot_display_gene_overview$cancer_type,levels=names(cancer_type.labs))




### plot ----

pdf(paste0(plot_dir,"oncoprint.clinrel_mutations_incl_prognostic.",plot_version,".pdf"),width=22,height=5)

p =ggplot(plot_display_gene_overview ,
       aes(x=patient_label,y=gene_name)) +
  geom_tile(aes(fill=alteration_display,color=clinrel_or_prognostic),size=0.5,width=0.8, height=0.8,alpha=1) + 
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(position="bottom") + xlab("") + ylab("") +
  geom_tile(data=filter(plot_display_gene_overview,gene_name=="annotation"),
            fill="white",size=1,width=0.7, height=0.7) +
  scale_color_manual(values = c(`TRUE`="black",`FALSE`="white"),na.value = "white") +
  scale_fill_manual(values=annotation_detail_colors,na.value="white") +
  facet_grid(cols=vars(cancer_type),scales="free",space = "free") 


p + facet_grid(cols=vars(cancer_type),scales="free",space = "free",labeller = labeller(cancer_type = cancer_type.labs)) 

p + facet_grid(cols=vars(cancer_type),scales="free",space = "free",labeller = labeller(cancer_type = cancer_type.labs)) + theme(legend.position="bottom")#, legend.direction = "vertical")

p + facet_rep_grid(cols=vars(cancer_type),scales="free",space = "free",repeat.tick.labels = "y") 



plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(gene_order))
plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_driver))


p =ggplot(plot_display_gene_overview ,
       aes(x=patient_label,y=gene_name)) +
  geom_tile(aes(fill=alteration_display,color=driver_gene),size=0.5,width=0.8, height=0.8,alpha=1) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(position="bottom") + xlab("") + ylab("") +
  geom_tile(data=filter(plot_display_gene_overview,gene_name=="annotation"),
            fill="white",size=1,width=0.7, height=0.7) +
  scale_color_manual(values = c(`TRUE`="black",`FALSE`="white"),na.value = "white") +
  scale_fill_manual(values=annotation_detail_colors,na.value="white") +
  facet_grid(cols=vars(cancer_type),scales="free",space = "free") 

p
p + facet_rep_grid(cols=vars(cancer_type),scales="free",space = "free",repeat.tick.labels = "y") 



dev.off()


## Add chrom alterations ----



chr_arm_cn_change_order = map_complex_sv_cn_change %>% select(chr_arm,call,chr_arm_cn_change) %>% arrange(chr_arm,call) %>% dplyr::mutate(chr_arm_cn_change=as.character(chr_arm_cn_change)) %>% select(chr_arm_cn_change) %>% unique() %>% flatten_chr()

map_complex_sv_cn_change$chr_arm_cn_change = factor(map_complex_sv_cn_change$chr_arm_cn_change,levels=chr_arm_cn_change_order)


#add to display gene
display_chr_arm_cn_change =  map_complex_sv_cn_change %>% 
  filter(chr_arm_cn_change %in% map_clinrel_chromalt$chr_arm_cn_change) %>% 
  mutate(gene_name=chr_arm_cn_change,alteration_display="complex",clinrel_or_prognostic=clinrel) %>% select(patient_label,cancer_type,gene_name,alteration_display,clinrel_or_prognostic)

display_chr_arm_cn_change_order = chr_arm_cn_change_order[chr_arm_cn_change_order %in% map_clinrel_chromalt$chr_arm_cn_change]


plot_display_gene_overview = rbind_no_colmatch(plot_display_gene_overview,display_chr_arm_cn_change)

plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(c(gene_order,display_chr_arm_cn_change_order)))

#order_patient_svs_cnt
plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_driver))





#add to display gene

display_unb_ctx =  unbalanced_ctx %>% filter(chr_arm_cn_change %in% map_clinrel_chromalt$chr_arm_cn_change) %>%  mutate(gene_name=chr_arm_cn_change,alteration_display="sv",clinrel_or_prognostic=clinrel) %>% select(patient_label,cancer_type,gene_name,alteration_display,clinrel_or_prognostic)

plot_display_gene_overview = rbind_no_colmatch(plot_display_gene_overview,display_unb_ctx)



all_cn_change = rbind(map_complex_sv_cn_change[,c("chr_arm","call","chr_arm_cn_change")],
                      unbalanced_ctx[,c("chr_arm","call","chr_arm_cn_change")])

chr_arm_cn_change_order = all_cn_change %>% select(chr_arm,call,chr_arm_cn_change) %>% dplyr::mutate(chr_arm_cn_change=as.character(chr_arm_cn_change)) %>% arrange(chr_arm,call) %>% select(chr_arm_cn_change) %>% unique() %>% flatten_chr()
display_chr_arm_cn_change_order = chr_arm_cn_change_order[chr_arm_cn_change_order %in% map_clinrel_chromalt$chr_arm_cn_change]




## plot unb ctx only [disabled] 

if(FALSE) {
pdf(paste0(plot_dir,"oncoprint.clinrel_mutations_incl_prognostic.chr_alt_svs.pdf"),width=22,height=8)

plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_effect))
plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(c(gene_order_incl_complex,display_chr_arm_cn_change_order)))


p =ggplot(plot_display_gene_overview ,
       aes(x=patient_label,y=gene_name)) +
  geom_tile(aes(fill=alteration_display,color=clinrel_or_prognostic),size=0.5,width=0.8, height=0.8,alpha=1) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(position="bottom") + xlab("") + ylab("") +
  geom_tile(data=filter(plot_display_gene_overview,gene_name=="annotation"),
            fill="white",size=1,width=0.7, height=0.7) +
  scale_color_manual(values = c(`TRUE`="black",`FALSE`="white"),na.value = "white") +
  scale_fill_manual(values=annotation_detail_colors,na.value="white") +
  facet_grid(cols=vars(cancer_type),scales="free",space = "free") 


p
p +  facet_rep_grid(cols=vars(cancer_type),scales="free",space = "free",repeat.tick.labels = "y") 


plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(c(gene_order,display_chr_arm_cn_change_order)))

#order_patient_svs_cnt
plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_driver))


p =ggplot(plot_display_gene_overview ,
       aes(x=patient_label,y=gene_name)) +
  geom_tile(aes(fill=alteration_display,color=clinrel_or_prognostic),size=0.5,width=0.8, height=0.8,alpha=1) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(position="bottom") + xlab("") + ylab("") +
  geom_tile(data=filter(plot_display_gene_overview,gene_name=="annotation"),
            fill="white",size=1,width=0.7, height=0.7) +
  scale_color_manual(values = c(`TRUE`="black",`FALSE`="white"),na.value = "white") +
  scale_fill_manual(values=annotation_detail_colors,na.value="white") +
  facet_grid(cols=vars(cancer_type),scales="free",space = "free") 

p
p +  facet_rep_grid(cols=vars(cancer_type),scales="free",space = "free",repeat.tick.labels = "y") 

dev.off()
}



#add numerical/aneuploidy


all_cn_change = rbind(map_complex_sv_cn_change[,c("chr_arm","call","chr_arm_cn_change")],
                      unbalanced_ctx[,c("chr_arm","call","chr_arm_cn_change")],
                      numerical_cn_chr_arms[,c("chr_arm","call","chr_arm_cn_change")])


chr_arm_cn_change_order = all_cn_change %>% select(chr_arm,call,chr_arm_cn_change) %>% dplyr::mutate(chr_arm_cn_change=as.character(chr_arm_cn_change)) %>% arrange(chr_arm,call) %>% select(chr_arm_cn_change) %>% unique() %>% flatten_chr()
display_chr_arm_cn_change_order = chr_arm_cn_change_order[chr_arm_cn_change_order %in% map_clinrel_chromalt$chr_arm_cn_change]


#add to display gene
display_numerical =  numerical_cn_chr_arms %>% filter(chr_arm_cn_change %in% clinrel_numerical_cn_chr_arms$chr_arm_cn_change) %>%  mutate(gene_name=chr_arm_cn_change,alteration_display="aneuploidy",clinrel_or_prognostic=clinrel) %>% select(patient_label,cancer_type,gene_name,alteration_display,clinrel_or_prognostic)

plot_display_gene_overview = rbind_no_colmatch(plot_display_gene_overview,display_numerical)

plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(c(gene_order,display_chr_arm_cn_change_order)))

#order_patient_svs_cnt
plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_driver))





## Final plot ----

pdf(paste0(plot_dir,"oncoprint.driver_genes.chrom_alt.",plot_version,".pdf"),width=22,height=8)
#pdf(paste0(plot_dir,"oncoprint.clinrel_mutations_incl_prognostic.chr_alt.",plot_version,".pdf"),width=22,height=8)
#pdf(paste0(plot_dir,"oncoprint.clinrel_mutations_incl_prognostic.chr_alt.",plot_version,".width18.pdf"),width=18,height=6)

plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_effect))
plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(c(gene_order_incl_complex,display_chr_arm_cn_change_order)))

p =ggplot(plot_display_gene_overview ,
       aes(x=patient_label,y=gene_name)) +
  geom_tile(aes(fill=alteration_display,color=clinrel_or_prognostic),size=0.5,width=0.8, height=0.8,alpha=1) + 
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(position="bottom") + xlab("") + ylab("") +
  geom_tile(data=filter(plot_display_gene_overview,gene_name=="annotation"),
            fill="white",size=1,width=0.7, height=0.7) +
  scale_color_manual(values = c(`TRUE`="black",`FALSE`="white"),na.value = "white") +
  scale_fill_manual(values=annotation_detail_colors_incl_aneuploidy,na.value="white") +
  facet_grid(cols=vars(cancer_type),scales="free",space = "free",labeller = labeller(cancer_type = cancer_type.labs))

p 

p  + theme(legend.position="bottom")#, legend.direction = "vertical")

p +  facet_rep_grid(cols=vars(cancer_type),scales="free",space = "free",repeat.tick.labels = "y",labeller = labeller(cancer_type = cancer_type.labs))  


plot_display_gene_overview$patient_label = factor(plot_display_gene_overview$patient_label, levels=(label_order_complex_driver))
plot_display_gene_overview$gene_name = factor(plot_display_gene_overview$gene_name,levels=rev(c(gene_order,display_chr_arm_cn_change_order)))


p =ggplot(plot_display_gene_overview ,
       aes(x=patient_label,y=gene_name)) +
  geom_tile(aes(fill=alteration_display,color=clinrel_or_prognostic),size=0.5,width=0.8, height=0.8,alpha=1) + 
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(position="bottom") + xlab("") + ylab("") +
  geom_tile(data=filter(plot_display_gene_overview,gene_name=="annotation"),
            fill="white",size=1,width=0.7, height=0.7) +
  scale_color_manual(values = c(`TRUE`="black",`FALSE`="white"),na.value = "white") +
  scale_fill_manual(values=annotation_detail_colors_incl_aneuploidy,na.value="white") +
  facet_grid(cols=vars(cancer_type),scales="free",space = "free",labeller = labeller(cancer_type = cancer_type.labs))

p 
p +  facet_rep_grid(cols=vars(cancer_type),scales="free",space = "free",repeat.tick.labels = "y",labeller = labeller(cancer_type = cancer_type.labs))  


dev.off()


#Export underlying data ----

#todo: narrow down export cols

export_cols = names(plot_display_gene_overview)
oncoprint_data = plot_display_gene_overview 
oncoprint_data_export = oncoprint_data %>% filter(gene_name!="annotation")
write.table(oncoprint_data_export[,export_cols],oncoprint_data_path,sep="\t",col.names=T,row.names=F,quote=F)

