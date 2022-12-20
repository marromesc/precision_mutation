# load libraries ----------------------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, reshape2, dplyr)
p_load(glue, janitor, params)
p_load(limma, broom)
p_load(forcats, effsize)
p_load(MultiAssayExperiment, SummarizedExperiment)
p_load(ComplexHeatmap, cowplot, ggsci, ggpubr, ggstatsplot)
p_load(LBSPR, magrittr,maftools, stargazer)
library(gridExtra);library(grid)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_indels/")
dir.create("02_filter_pindel_files_v1")
odir = "02_filter_pindel_files_v1"
setwd("~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_indels/02_filter_pindel_files_v1/")

data_folder = "~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/analysis_indels/01_read_pindel_files_v1/"

source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/mut2maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/to_vcf.maf.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_sheets.R')
source('~/sea_rsrch3_home/rsrch2_backup/pers_docs/sahil_scripts/gm_mutect.R')
install.packages('~/sea_rsrch3_home/rsrch2_backup/projects/DCIS_exome/scripts_run2/github_wranglr/', repos = NULL, type="source")
library(wranglr)

# read files ----------
df_indels = read_rds(path = paste0(data_folder, "df_indels5.rds"))
dim(df_indels) #[1] 431034    117

df_indels[df_indels == "."] <- NA
table(is.na(df_indels))

# FALSE     TRUE 
# 24762897 25668081 

df_indels$variant_classification = recode_variant_classification(df_indels$exonicfunc.knowngene)
df_indels = df_indels %>% mutate(key = paste0(chrom, ":", start, "-", end, ":", type, "-", length)) %>% dplyr::mutate(aaannotation1 = gsub("[p[:punct:]]", "", aaannotation, ignore.case = FALSE),
                                                                                                                            aaannotation2 = gsub("([A-Z]?)([0-9]{,6})(.*)", "\\1\\2", aaannotation1),
                                                                                                                            hotkey = paste0(gene.knowngene, "-", aaannotation2),
                                                                                                                            normal_af = n_alt_count/(n_ref_count + n_alt_count),
                                                                                                                            t_depth = t_ref_count + t_alt_count,
                                                                                                                            n_depth = n_ref_count + n_alt_count)

df_indels$exac_all = as.numeric(df_indels$exac_all)
df_indels$esp6500siv2_all = as.numeric(df_indels$esp6500siv2_all)
df_indels$x1kg2015aug_max = as.numeric(df_indels$x1kg2015aug_max)
df_indels$t_depth = as.numeric(df_indels$t_depth)
df_indels$n_depth = as.numeric(df_indels$n_depth)
df_indels$normal_af = as.numeric(df_indels$normal_af)
df_indels$tumor_f = as.numeric(df_indels$tumor_f)
df_indels$t_vaf = as.numeric(df_indels$t_vaf)
df_indels$n_vaf = as.numeric(df_indels$n_vaf)

table(is.na(df_indels$t_vaf))
table(is.na(df_indels$n_vaf))
table(is.na(df_indels$t_depth))
table(is.na(df_indels$n_depth))

# ** remove samples in exac databases  --------
df_indels = df_indels %>% dplyr::filter(t_vaf >= 0.02, esp6500siv2_all < 0.01 | is.na(esp6500siv2_all), 
                                                          exac_all < 0.01 | is.na(exac_all), 
                                                          x1kg2015aug_max < 0.01 | is.na(x1kg2015aug_max))

# ** remove samples with low coverage  --------
df_indels = df_indels %>% dplyr::filter(t_depth >= 15, n_depth >= 10, n_vaf <= 0.01)
dim(df_indels) # 73295   125

# ** keep exonic SNV's --------
df_indels = df_indels %>% dplyr::filter(func.knowngene %in% c("exonic"))
dim(df_indels) # 60406   125

df_indels2 = unique(df_indels)
dim(df_indels2)


write.csv(df_indels2, "df_indels2.csv")
df_indels2 = read.csv("df_indels2.csv")

# ** remove TD's  --------
df_indels2 = df_indels2 %>% dplyr::filter(!type %in% c("TD"))
dim(df_indels2) # 6417   125

df_indels2 = unique(df_indels2)
dim(df_indels2)

write.csv(df_indels2, "df_indels2_withoutTDs.csv")

# count the indels per sample -----
df_indels_notds = df_indels2 %>% dplyr::group_by(sample_id) %>% mutate(mut_cnt = n())
write.csv(df_indels_notds, "df_indels_notds.csv")

df_indels_cnt = df_indels_notds %>% dplyr::distinct(sample_id, mut_cnt)
p1 = df_indels_cnt
pdf("pindel_count_per_sample.pdf", height = 42)
p<-tableGrob(p1)
grid.arrange(p)
dev.off()

df_indels_notds$matching_key = paste0(df_indels_notds$sample_id, "_", df_indels_notds$key, "_", df_indels_notds$hotkey)
write.csv(df_indels_notds, "df_indels_notds.csv")
df_indels_notds = read.csv("df_indels_notds.csv")
df_indels_notds = df_indels_notds[,-1]


# ** filter for strand bias ------
df_indels_notds1 = df_indels_notds %>% dplyr::filter(supportup >= 2, supportdown >= 2)

# annotate shared and private ----
df_merged10 = df_indels_notds1 %>% group_by(key, patient) %>% mutate(n_index = n(), shared = case_when(n_index >= 2 ~ "shared", n_index == 1 ~ "private")) %>% ungroup()
write.csv(df_merged10, "df_merged10.csv")
dim(df_merged10)
#[1] 3202  130

# count the indels per sample -----
df_merged10cnt = df_merged10 %>% dplyr::group_by(sample_id) %>% mutate(mut_cnt = n())

df_merged10cnt = df_merged10cnt %>% dplyr::distinct(sample_id, mut_cnt)

p2 = df_merged10cnt
pdf("pindel_count_per_sample_post_strandfilter2.pdf", height = 42)
p<-tableGrob(p2)
grid.arrange(p)
dev.off()


# count the shared indels per patient -----
df_merged11cnt = df_merged10 %>% dplyr::group_by(patient, shared) %>% mutate(mut_cnt = n())

df_merged11cnt = df_merged11cnt %>% dplyr::distinct(patient, shared, mut_cnt)

p3 = df_merged11cnt
pdf("pindel_count_per_patient_post_strandfilter2.pdf", height = 36)
p<-tableGrob(p3)
grid.arrange(p)
dev.off()
